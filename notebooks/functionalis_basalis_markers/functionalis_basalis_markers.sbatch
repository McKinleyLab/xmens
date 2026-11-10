#!/bin/bash
# Runs functionalis_basalis_markers.ipynb with papermill. Submit it from the folder that should receive the outputs:
#
#   HERE=/n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks/functionalis_basalis_markers
#   mkdir -p /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_08_functionalis_basalis_markers
#   cd /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_08_functionalis_basalis_markers
#   sbatch $HERE/functionalis_basalis_markers.sbatch
#
# It writes there the executed notebook, strict_pairs.csv, the two shortlists, figures/ and the Slurm log, and stops
# before writing anything if figures/ is already there. A run takes about 6 minutes and 4.2 GB.
#SBATCH -J fb_markers
#SBATCH -p eddy
#SBATCH --account=eddy_lab
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o functionalis_basalis_markers-%j.out
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
set -euo pipefail

NH=/n/lab_storage/eddy_lab/users/nhilgert
HERE=$NH/xmens_scripts/notebooks/functionalis_basalis_markers
CACHE=$NH/.cache_scratch

# Nothing may be written in the home folder: caches and Jupyter and IPython state go to lab storage.
export HOME=$NH/tmp
export MPLCONFIGDIR=$CACHE/mpl NUMBA_CACHE_DIR=$CACHE/numba XDG_CACHE_HOME=$CACHE PYTHONPYCACHEPREFIX=$CACHE/pycache
export IPYTHONDIR=$CACHE/ipython JUPYTER_DATA_DIR=$CACHE/jupyter_data JUPYTER_CONFIG_DIR=$CACHE/jupyter_config
export JUPYTER_RUNTIME_DIR=$NH/tmp/jupyter_runtime_functionalis_basalis_markers
export HDF5_USE_FILE_LOCKING=FALSE

# One BLAS thread: OpenBLAS otherwise starts a thread for every core of the node, not of the job.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$JUPYTER_RUNTIME_DIR" "$PYTHONPYCACHEPREFIX"

# --cwd: the notebook's relative paths start from its own folder.
"$NH/.uv/envs/rsc_uv/bin/python" -m papermill \
    --kernel python3 \
    --cwd "$HERE" \
    --log-output \
    "$HERE/functionalis_basalis_markers.ipynb" \
    functionalis_basalis_markers.executed.ipynb \
    -p OUT_DIR "$PWD"
