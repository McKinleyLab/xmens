# Step 2: Monocle 2 (monocle 2.34.0) from the raw counts to the z matrix of every pseudotime heatmap of one species.
#   1. newCellDataSet() on the counts of every gene for the nuclei with a pseudotime (step 1), negbinomial.size()
#   2. estimateSizeFactors() and estimateDispersions(), Monocle's defaults, on all of those nuclei and genes
#   3-4. plot_pseudotime_heatmap()'s own lines, copied verbatim (monocle_z below): genSmoothCurves(), a negative-
#      binomial fit of the size-factor-normalized counts on sm.ns(Pseudotime, df = 3) with the dispersion from the
#      trend, at 100 pseudotime points; then log10(curve + 1) and z per gene. Its last step, the clip to +/-3, is
#      applied when drawing (step 3), as it is in Monocle, right before the colours
# Mouse also gives Fig. 4H: the Yang et al. (2023) panel on the induced nuclei, pseudotime re-ranked within them (0 to
# 1), with the size factors and dispersion trend of all nuclei; and plot_pseudotime_heatmap()'s own drawing of it.
# Writes into OUT_DIR:
#   z_<species>.tsv                 gene, fit, then z at the 100 points, for every MI-heatmap and Fig. 5F gene; fit is
#                                   "trend", or "fallback" for a gene whose fit failed and Monocle refitted with
#                                   negbinomial()
#   dispersion_<species>.tsv.gz     dispersionTable(): per gene mean, fitted and empirical dispersion
#   dispersion_<species>_trend.tsv  the trend's coefficients, asymptDisp a and extraPois b
#   sessionInfo_<species>.txt
#   mouse only: fig4h_z.tsv (rows named as published), fig4h_monocle_plot_pseudotime_heatmap.pdf
# usage: Rscript monocle_z.R SPECIES INPUTS_DIR OUT_DIR
suppressPackageStartupMessages({
    library(monocle)
    library(Matrix)
})
# estimateDispersions() calls dplyr functions that dplyr 1.2 made defunct; step 0 installs dplyr 1.1.4 beside monocle.
if (packageVersion("dplyr") >= "1.2.0") {
    stop("dplyr ", packageVersion("dplyr"), " is loaded; Monocle 2 needs dplyr 1.1.4 (step 0), with its library first ",
         "on the search path (R_LIBS=$MONOCLE_LIB)")
}

args <- commandArgs(trailingOnly = TRUE)
species <- args[1]
inputs <- args[2]
out_dir <- args[3]
stopifnot(species %in% c("mouse", "human"))
if (file.exists(file.path(out_dir, sprintf("z_%s.tsv", species)))) {
    stop(file.path(out_dir, sprintf("z_%s.tsv", species)), " exists; not overwriting it")
}
# Files are written into a temporary folder and moved into OUT_DIR when all are complete, z_<species>.tsv last.
partial <- file.path(out_dir, paste0(".partial_", species))
dir.create(partial, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
say <- function(...) cat(sprintf("[%4.0f s] ", as.numeric(difftime(Sys.time(), started, units = "secs"))), ..., "\n",
                         sep = "")
write_tsv <- function(x, name) write.table(x, file.path(partial, name), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 1-2. the CellDataSet, size factors and dispersions --------------------------------------------------------------
counts <- as(readMM(file.path(inputs, species, "counts.mtx")), "CsparseMatrix")       # genes x nuclei
genes <- readLines(file.path(inputs, species, "genes.txt"))
cells <- read.delim(file.path(inputs, sprintf("cells_%s.tsv.gz", species)))
dimnames(counts) <- list(genes, cells$cell)
cds <- newCellDataSet(counts,
                      phenoData = new("AnnotatedDataFrame", data = data.frame(Pseudotime = cells$pt_rank_norm,
                                                                               timepoint = cells$timepoint,
                                                                               row.names = cells$cell)),
                      featureData = new("AnnotatedDataFrame", data = data.frame(gene_short_name = genes,
                                                                                 row.names = genes)),
                      expressionFamily = negbinomial.size())
say(species, ": ", ncol(cds), " nuclei x ", nrow(cds), " genes")
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
coefficients <- attr(cds@dispFitInfo[["blind"]]$disp_func, "coefficients")
say("dispersion trend over ", nrow(cds@dispFitInfo[["blind"]]$disp_table), " genes: asymptDisp ",
    signif(coefficients[1], 6), ", extraPois ", signif(coefficients[2], 6))
write_tsv(data.frame(coefficient = names(coefficients), value = unname(coefficients)),
          sprintf("dispersion_%s_trend.tsv", species))
gz <- gzfile(file.path(partial, sprintf("dispersion_%s.tsv.gz", species)), "w")
write.table(dispersionTable(cds), gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)

# ---- 3-4. curves and z, as plot_pseudotime_heatmap() makes them ------------------------------------------------------
monocle_z <- function(cds_subset) {
    # Between the markers: plot_pseudotime_heatmap()'s lines (monocle 2.34.0), verbatim, from its call of
    # genSmoothCurves() to z, with its defaults (trend_formula, cores, and norm_method = "log", of whose if/else only
    # the log10 line runs); its row clustering and drawing are left out. Step 0's check finds each line in the
    # installed function. Monocle prints the error of a gene whose fit fails before refitting it with negbinomial():
    # the output is captured, and if anything was printed, each gene is refitted alone to name those genes.
    trend_formula <- "~sm.ns(Pseudotime, df=3)"
    cores <- 1
    printed <- capture.output({
        # >>> plot_pseudotime_heatmap()
        pseudocount <- 1
        newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), max(pData(cds_subset)$Pseudotime), length.out = 100))
        m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula, relative_expr = T, new_data = newdata)
        m = m[!apply(m, 1, sum) == 0, ]
        m = log10(m + pseudocount)
        m = m[!apply(m, 1, sd) == 0, ]
        m = Matrix::t(scale(Matrix::t(m), center = TRUE))
        m = m[is.na(row.names(m)) == FALSE, ]
        m[is.nan(m)] = 0
        # <<< its next two lines, m[m > scale_max] = scale_max and m[m < scale_min] = scale_min (+/-3), are applied
        #     when drawing (step 3)
    })
    fit <- rep("trend", nrow(m))
    if (length(printed) > 0) {
        fit <- vapply(rownames(cds_subset), function(g) {
            alone <- capture.output(genSmoothCurves(cds_subset[g, ], cores = cores, trend_formula = trend_formula,
                                                    relative_expr = T, new_data = newdata))
            if (length(alone) > 0) "fallback" else "trend"
        }, character(1))
    }
    if (nrow(m) != nrow(cds_subset)) stop("genes without a curve: ", paste(setdiff(rownames(cds_subset), rownames(m)),
                                                                           collapse = ", "))
    data.frame(gene = rownames(m), fit = fit, as.matrix(m), check.names = FALSE)
}

heatmap_genes <- unique(c(read.csv(file.path(inputs, sprintf("mi_top_%s.csv", species)))$gene,
                          read.csv(file.path(inputs, "fig5f_genes.csv"))[[species]]))
z <- monocle_z(cds[heatmap_genes, ])
say(nrow(z), " genes on ", ncol(cds), " nuclei; fallback fits: ", sum(z$fit != "trend"))

if (species == "mouse") {
    yang <- read.delim(file.path(inputs, "yang_panel.tsv"))
    induced <- pData(cds)$timepoint != "uninduced"
    cds_yang <- cds[yang$symbol, induced]
    pData(cds_yang)$Pseudotime <- (rank(pData(cds_yang)$Pseudotime, ties.method = "first") - 1) / (sum(induced) - 1)
    fData(cds_yang)$gene_short_name <- yang$gene                     # published names
    z_yang <- monocle_z(cds_yang)
    z_yang$gene <- yang$gene
    say(nrow(z_yang), " Yang panel genes on ", ncol(cds_yang), " induced nuclei; fallback fits: ",
        sum(z_yang$fit != "trend"))
    write_tsv(z_yang, "fig4h_z.tsv")
    pdf(file.path(partial, "fig4h_monocle_plot_pseudotime_heatmap.pdf"), width = 3, height = 4)
    plot_pseudotime_heatmap(cds_yang, cluster_rows = FALSE, show_rownames = TRUE)   # Monocle's own drawing
    invisible(dev.off())
}

writeLines(capture.output(sessionInfo()), file.path(partial, sprintf("sessionInfo_%s.txt", species)))
done <- sprintf("z_%s.tsv", species)
write_tsv(z, done)
for (f in c(setdiff(list.files(partial), done), done)) file.rename(file.path(partial, f), file.path(out_dir, f))
unlink(partial, recursive = TRUE)
say("wrote ", out_dir)
