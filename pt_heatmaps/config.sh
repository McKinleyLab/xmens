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

# The final objects (read only): fibroblast and decidual nuclei, raw counts in layers['counts'], and the rank-normalized
# pseudotime in obs['pseudotime'] (NaN for the 907 apoptotic mouse decidual nuclei). netscratch is purged after 90 days:
# point OBJECTS at a permanent copy of the same files once there is one.
OBJECTS=/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects
H5AD_MOUSE=$OBJECTS/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad
H5AD_HUMAN=$OBJECTS/Cevrim_Human_Only_Fib_Dec_Cells.h5ad

# The MI of record and the one-to-one ortholog pairs, which give the genes of the big MI heatmaps, and the 28 genes of
# Fig. 5F in drawing order. xmens_mi makes the same files ($OUT of xmens_mi: mi_mouse.csv, mi_human.csv,
# overlap/ortholog_pairs_one2one.csv and fig5f/conserved_profiles_selected_genes.csv).
RESULTS=$NH/xmens_final/results/mi_dpt
MI_MOUSE=$RESULTS/pt_harmonized_mouse_nall/mi_all.csv
MI_HUMAN=$RESULTS/pt_harmonized_human_nall/mi_all.csv
PAIRS=$RESULTS/overlap_of_record/ortholog_pairs_one2one.csv
FIG5F_GENES=$RESULTS/overlap_of_record/crossspecies_heatmap_selection/conserved_profiles_selected_genes.csv

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
