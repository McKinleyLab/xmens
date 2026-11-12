# Paths and settings shared by every step; each step script sources this file.
#
# OUT, the run folder, has no default, so a run never writes into an earlier one by accident. Set it once per run:
#   export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_pt_heatmaps

OUT=${OUT:?export OUT, the run folder, before running a step (see README.md)}
NH=/n/lab_storage/eddy_lab/users/nhilgert

# Python (rsc_uv: anndata, numpy, pandas, scipy, matplotlib 3.10) for steps 1 and 3; R 4.4.1 for step 2.
PY=$NH/.uv/envs/rsc_uv/bin/python
R_MODULE=R/4.4.1-fasrc01

# Monocle 2 (Bioconductor 3.20, monocle 2.34.0), the dependencies the lab's R library lacks, and dplyr 1.1.4 (Monocle's
# estimateDispersions() fails with the lab's dplyr 1.2.1), in their own library (step 0); step 2 puts it first on R's
# search path. The lab's main library, $R_LIBS_USER, supplies everything else and is not changed.
MONOCLE_LIB=$NH/R/library-4.4.1-monocle2

# The final objects (read only): the fibroblast/decidual objects of 11 Nov 2026 with the pseudotime of 10 Nov in
# obs['pseudotime'] (NaN for the apoptotic mouse decidual nuclei) and raw counts in layers[$COUNTS_LAYER]. (The
# deposited objects of 2 Nov, drawn in the record, are $NS/final_objects/Cevrim_*_Only_Fib_Dec_Cells.h5ad with the
# counts in layers['counts'].)
OBJECTS=/n/lab_storage/eddy_lab/Lab/mckinley/data
H5AD_MOUSE=$OBJECTS/20261111_mouse_only_eFib_dec_final_nick.h5ad
H5AD_HUMAN=$OBJECTS/20261111_Human_only_eFib_dec_final_nick.h5ad
COUNTS_LAYER=raw_counts
# The MI run of November 2026 (xmens_mi on the objects above) and its one-to-one ortholog pairs, which give the genes
# of the big MI heatmaps, and its 28 genes of Fig. 5F in drawing order. (The MI of record, drawn in the record:
# $NH/xmens_final/results/mi_dpt/pt_harmonized_{mouse,human}_nall/mi_all.csv and overlap_of_record/.)
RESULTS=/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_11_xmens_mi_november
MI_MOUSE=$RESULTS/mi_mouse.csv
MI_HUMAN=$RESULTS/mi_human.csv
PAIRS=$RESULTS/overlap/ortholog_pairs_one2one.csv
FIG5F_GENES=$RESULTS/fig5f/conserved_profiles_selected_genes.csv
# Caches, temporary files and R packages go to lab storage, never to the home directory.
CACHE=$NH/.cache_scratch
export HDF5_USE_FILE_LOCKING=FALSE
export MPLCONFIGDIR=$CACHE/mpl
export NUMBA_CACHE_DIR=$CACHE/numba
export XDG_CACHE_HOME=$CACHE/xdg
export IPYTHONDIR=$CACHE/ipython                # a notebook run's kernel history
export JUPYTER_RUNTIME_DIR=$CACHE/jupyter_runtime  # and its connection files
export R_LIBS_USER=$NH/R/library-4.4.1
export R_MAKEVARS_USER=$NH/R/Makevars
export TMPDIR=$NH/tmp

# Prints "found" or "missing", for --check.
exists() {
    if [ -e "$1" ]; then
        echo found
    else
        echo missing
    fi
}
