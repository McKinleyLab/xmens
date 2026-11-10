# Paths and environment for the DESeq2 notebooks; every sbatch script in this folder sources it.
# $HERE (this folder) is set by the sbatch script. Export OUT, the run folder, first (see README.md).
OUT=${OUT:?export OUT, the run folder, before running a notebook (see README.md)}

NH=/n/lab_storage/eddy_lab/users/nhilgert
NS=/n/netscratch/eddy_lab/Everyone/nhilgert
REPO=$NH/xmens_final

# Inputs of record: the objects the published fits were made on.
XMENS_H5AD=$REPO/data/interim/final_mouse.h5ad
HUMAN_H5AD=$REPO/data/deliverables/2026_09_15/Cevrim_Human_Only_Fib_Dec_Cells.h5ad
ULRICH_H5AD=$NS/ulrich2024/ulrich2024_20donor_stromal.h5ad

# Small inputs and the tables of record, copied into this folder.
RESOURCES=$HERE/resources
RECORD=$HERE/record

PY=$NH/.uv/envs/rsc_uv/bin/python
PAPERMILL=$NH/.uv/envs/rsc_uv/bin/papermill
export JUPYTER_PATH=$HERE/jupyter

# Caches and the Jupyter and IPython runtime folders go to lab storage: nothing may be written in the home folder.
CACHE=$NH/.cache_scratch
export HDF5_USE_FILE_LOCKING=FALSE
export MPLCONFIGDIR=$CACHE/mpl NUMBA_CACHE_DIR=$CACHE/numba XDG_CACHE_HOME=$CACHE/xdg
export IPYTHONDIR=$CACHE/ipython JUPYTER_RUNTIME_DIR=$CACHE/jupyter_runtime
export JUPYTER_DATA_DIR=$CACHE/jupyter_data JUPYTER_CONFIG_DIR=$CACHE/jupyter_config
# The rsc_uv interpreter lives under the home folder (~/.local/share/uv/python), so bytecode it compiles is kept here.
export PYTHONPYCACHEPREFIX=$CACHE/pycache
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$XDG_CACHE_HOME" "$IPYTHONDIR" "$JUPYTER_RUNTIME_DIR" "$PYTHONPYCACHEPREFIX"

# PyDESeq2 parallelizes over genes with joblib (N_CPUS); BLAS threads on top of that only compete for the same cores.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Runs a notebook with papermill into a .partial copy and renames it when the run completes, so a failed run never
# looks like a finished one.
run_notebook() {
    local notebook=$1 kernel=$2 executed=$3
    shift 3
    local partial=${executed%.ipynb}.partial.ipynb
    mkdir -p "$(dirname "$executed")"
    "$PAPERMILL" "$notebook" "$partial" -k "$kernel" --cwd "$(dirname "$notebook")" --log-output "$@"
    mv "$partial" "$executed"
}

exists() {
    if [ -e "$1" ]; then
        echo found
    else
        echo missing
    fi
}
