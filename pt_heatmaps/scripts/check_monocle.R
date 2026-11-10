# Step 0 check, run as step 2 runs R (the library from install_monocle.R first on the search path, through R_LIBS):
# monocle 2.34.0 loads with dplyr 1.1.4; estimateSizeFactors() and estimateDispersions() run on 200 simulated genes
# x 60 cells; and every line monocle_z.R copies from plot_pseudotime_heatmap() is a line of the installed function.
# Writes <library>/monocle_install_sessionInfo.txt.
# usage: R_LIBS=LIBRARY Rscript check_monocle.R LIBRARY
lib <- commandArgs(trailingOnly = TRUE)[1]
suppressPackageStartupMessages(library(monocle))
for (p in c("monocle", "dplyr", "VGAM", "Biobase")) {
    cat(p, as.character(packageVersion(p)), "loads from", dirname(find.package(p)), "\n")
}
stopifnot(packageVersion("monocle") == "2.34.0", packageVersion("dplyr") == "1.1.4")

set.seed(1)
counts <- matrix(rnbinom(200 * 60, mu = rep(2^seq(-2, 6, length.out = 200), 60), size = 2), 200, 60,
                 dimnames = list(paste0("gene", 1:200), paste0("cell", 1:60)))
test <- newCellDataSet(counts, featureData = new("AnnotatedDataFrame",
                                                 data = data.frame(gene_short_name = rownames(counts),
                                                                   row.names = rownames(counts))),
                       expressionFamily = negbinomial.size())
test <- estimateDispersions(estimateSizeFactors(test))
cat("estimateDispersions() ran on simulated data: dispersion trend over", nrow(dispersionTable(test)), "genes\n")

# The lines monocle_z.R copies from plot_pseudotime_heatmap() (between its >>> and <<< markers) must each be a line
# of the installed function.
ours <- readLines(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "monocle_z.R"))
ours <- trimws(ours[(grep(">>> plot_pseudotime_heatmap", ours) + 1):(grep("<<<", ours) - 1)])
theirs <- trimws(deparse(body(plot_pseudotime_heatmap), width.cutoff = 500))
missing <- setdiff(ours, theirs)
if (length(missing) > 0) stop("not in plot_pseudotime_heatmap(): ", paste(missing, collapse = " | "))
cat("monocle_z.R:", length(ours), "lines copied from plot_pseudotime_heatmap(), each found in the installed function\n")

writeLines(capture.output(sessionInfo()), file.path(lib, "monocle_install_sessionInfo.txt"))
