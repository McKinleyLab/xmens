#!/bin/bash
# Runs stro7_cosine.ipynb with papermill. Submit it from the folder that should receive the outputs, for example
#
#   mkdir -p /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_08_stro7_cosine
#   cd /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_08_stro7_cosine
#   sbatch /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/stro7_cosine/stro7_cosine.sbatch
#
# It writes there the executed notebook, stro7_cosine_human.csv, stro7_cosine_mouse.csv and the Slurm log.
#SBATCH -J stro7_cosine
#SBATCH -p eddy
#SBATCH --account=eddy_lab
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o stro7_cosine-%j.out
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
set -euo pipefail

NH=/n/lab_storage/eddy_lab/users/nhilgert
CACHE=$NH/.cache_scratch

# Nothing may be written in the home folder: caches and Jupyter and IPython state go to lab storage.
export HOME=$NH/tmp
export MPLCONFIGDIR=$CACHE/mpl NUMBA_CACHE_DIR=$CACHE/numba XDG_CACHE_HOME=$CACHE PYTHONPYCACHEPREFIX=$CACHE/pycache
export IPYTHONDIR=$CACHE/ipython JUPYTER_DATA_DIR=$CACHE/jupyter_data JUPYTER_CONFIG_DIR=$CACHE/jupyter_config
export JUPYTER_RUNTIME_DIR=$NH/tmp/jupyter_runtime_stro7_cosine
export HDF5_USE_FILE_LOCKING=FALSE

# One BLAS thread: OpenBLAS otherwise starts a thread for every core of the node, not of the job.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$JUPYTER_RUNTIME_DIR" "$PYTHONPYCACHEPREFIX"

"$NH/.uv/envs/rsc_uv/bin/python" -m papermill \
    --kernel python3 \
    --log-output \
    "$NH/xmens_scripts/stro7_cosine/stro7_cosine.ipynb" \
    stro7_cosine.executed.ipynb \
    -p OUT_DIR "$PWD"
