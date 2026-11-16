# Step 0: install Monocle 2 from Bioconductor 3.20 (the release for R 4.4; monocle 2.34.0) into the library given as the
# first argument, with only the dependencies missing from every library on .libPaths(). Packages that are already
# installed elsewhere are used as they are, and nothing is updated, so the lab's main library is not changed.
# One exception: dplyr 1.1.4 also goes into this library. estimateDispersions() calls dplyr's group_by_() and select_(),
# which dplyr 1.2.0 made defunct (they stop with an error; the lab's library has dplyr 1.2.1); in dplyr 1.1.4 they still
# run, with a deprecation warning. Scripts that put this library first on the search path (R_LIBS) load dplyr 1.1.4.
# usage: Rscript install_monocle.R LIBRARY
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
cat("libraries:", .libPaths(), sep = "\n  ")
options(Ncpus = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

# Nothing is loaded here (loading monocle would load the lab's dplyr 1.2.1); scripts/check_monocle.R then loads monocle
# in a new R session, as step 2 does.
if (nzchar(system.file(package = "monocle"))) {
    cat("monocle", as.character(packageVersion("monocle")), "is already installed in",
        dirname(system.file(package = "monocle")), "\n")
} else {
    BiocManager::install("monocle", version = "3.20", lib = lib, update = FALSE, ask = FALSE)
}
if (dir.exists(file.path(lib, "dplyr")) && packageVersion("dplyr", lib.loc = lib) == "1.1.4") {
    cat("dplyr 1.1.4 is already installed in", lib, "\n")
} else {
    install.packages("https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.1.4.tar.gz", repos = NULL,
                     type = "source", lib = lib)
}
