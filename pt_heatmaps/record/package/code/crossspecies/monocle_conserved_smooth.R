# Monocle 2 plot_pseudotime_heatmap() smoother, one species per run (pseudotime heatmaps of the X-Mens revision).
# Fitting code as in Monocle 2 (parametricDispersionFit with Cook's-distance refit; per gene vglm of size-factor-
# normalized counts on sm.ns(Pseudotime, df = 3) with negbinomial.size(size = 1 / dispersion trend), epsilon = 1e-1;
# predict(type = "response") at 100 evenly spaced pseudotime points). Inputs:
#   * PERCELL: one row per nucleus with total_counts, pt_rank_norm (rank-normalized pseudotime) and cnt_<gene> columns;
#     the genes smoothed are the cnt_* columns; size factors = total_counts / geometric mean over all rows
#   * DISP_TABLE: per-gene (mu, disp) for every gene detected in > 1 nucleus (Monocle estimateDispersions behaviour);
#     the dispersion-mean trend is fitted across all of them
# Outputs the raw response curves (OUT_TSV: gene, fit, 100 columns) and the grid (<OUT>_grid.tsv); log10(x + 1),
# z-scoring and clipping are done at render time.
# usage: Rscript monocle_conserved_smooth.R PERCELL_TSV_GZ DISP_TABLE_TSV_GZ OUT_TSV NCORES
suppressPackageStartupMessages({library(VGAM); library(splines); library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
percell <- args[1]; disp_tsv <- args[2]; out_tsv <- args[3]; ncores <- as.integer(args[4])

pc <- read.delim(gzfile(percell), check.names = FALSE)
cc <- grep("^cnt_", colnames(pc), value = TRUE); genes <- sub("^cnt_", "", cc)
cnt <- as.matrix(pc[, cc]); colnames(cnt) <- genes
sf_all <- pc$total_counts / exp(mean(log(pc$total_counts)))

disp_table <- read.delim(gzfile(disp_tsv))
rownames(disp_table) <- make.unique(disp_table$gene)
disp_table <- disp_table[is.finite(disp_table$mu) & disp_table$mu > 0, c("mu", "disp")]
disp_table$disp[!is.finite(disp_table$disp) | disp_table$disp < 0] <- 0

parametricDispersionFit <- function(disp_table, initial_coefs = c(1e-6, 1)) {
  coefs <- initial_coefs; iter <- 0
  while (TRUE) {
    residuals <- disp_table$disp / (coefs[1] + coefs[2] / disp_table$mu)
    good <- disp_table[which((residuals > initial_coefs[1]) & (residuals < 10000)), ]
    fit <- suppressWarnings(glm(disp ~ I(1/mu), data = good, family = Gamma(link = "identity"), start = coefs))
    oldcoefs <- coefs; coefs <- coefficients(fit)
    if (coefs[1] < initial_coefs[1]) coefs[1] <- initial_coefs[1]
    if (coefs[2] < 0) stop("Parametric dispersion fit failed")
    if (sum(log(coefs / oldcoefs)^2) < initial_coefs[1]) break
    iter <- iter + 1; if (iter > 10) { warning("Dispersion fit did not converge."); break }
  }
  list(fit, coefs)
}
res <- parametricDispersionFit(disp_table)
CD <- cooks.distance(res[[1]]); outl <- union(names(CD[CD > 4 / nrow(disp_table)]), setdiff(rownames(disp_table), names(CD)))
res <- parametricDispersionFit(disp_table[!rownames(disp_table) %in% outl, ])
coefs <- res[[2]]
disp_func <- function(q) coefs[1] + coefs[2] / q
cat(sprintf("dispersion trend over %d genes: asymptDisp=%.4g extraPois=%.4g (outliers removed: %d)\n",
            nrow(disp_table), coefs[1], coefs[2], length(outl)))

Pseudotime <- pc$pt_rank_norm
newdata <- data.frame(Pseudotime = seq(min(Pseudotime), max(Pseudotime), length.out = 100))
sf3 <- sf_all
smooth_gene <- function(x) {
  f_expression <- round(x / sf3)
  fam <- negbinomial.size(size = 1 / disp_func(mean(round(x))))
  used <- "trend"
  fit <- tryCatch(suppressWarnings(vglm(f_expression ~ sm.ns(Pseudotime, df = 3), family = fam, epsilon = 1e-1)),
                  error = function(e) { used <<- "fallback_negbinomial"
                    suppressWarnings(vglm(f_expression ~ sm.ns(Pseudotime, df = 3), family = negbinomial(), epsilon = 1e-1, checkwz = TRUE)) })
  list(curve = as.numeric(predict(fit, newdata = newdata, type = "response")), used = used)
}
t0 <- Sys.time()
fits <- mclapply(genes, function(g) smooth_gene(cnt[, g]), mc.cores = ncores)
bad <- vapply(fits, function(f) inherits(f, "try-error") || is.null(f$curve), logical(1))
if (any(bad)) stop(sprintf("%d genes failed: %s", sum(bad), paste(genes[bad], collapse = ",")))
m <- t(vapply(fits, function(f) f$curve, numeric(100)))
used <- vapply(fits, function(f) f$used, character(1))
cat(sprintf("fit %d genes x %d cells in %.1f s; fallback used for %d: %s\n", length(genes), nrow(pc),
            as.numeric(difftime(Sys.time(), t0, units = "secs")), sum(used != "trend"), paste(genes[used != "trend"], collapse = ",")))
write.table(data.frame(gene = genes, fit = used, m, check.names = FALSE), out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(newdata, sub("\\.tsv$", "_grid.tsv", out_tsv), sep = "\t", quote = FALSE, row.names = FALSE)
cat("done\n")
