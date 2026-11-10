#!/bin/bash
#SBATCH -J pt_mouse
#SBATCH -p eddy
#SBATCH --account=eddy_lab
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 03:00:00
#SBATCH --constraint=icelake
#SBATCH -o /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks/pseudotime_runs/pt_mouse-%j.out
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
set -euo pipefail
C=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch
export HOME=/n/lab_storage/eddy_lab/users/nhilgert/tmp MPLCONFIGDIR=$C/mpl NUMBA_CACHE_DIR=$C/numba IPYTHONDIR=$C/ipython XDG_CACHE_HOME=$C \
       JUPYTER_DATA_DIR=$C/jupyter_data JUPYTER_CONFIG_DIR=$C/jupyter_config JUPYTER_RUNTIME_DIR=/n/lab_storage/eddy_lab/users/nhilgert/tmp/jupyter_runtime_pt JAX_PLATFORMS=cpu
mkdir -p $JUPYTER_RUNTIME_DIR
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks
hostname; lscpu | grep 'Model name'; nproc
/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python -m papermill --log-output --kernel python3 pseudotime_mouse.ipynb /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks/pseudotime_runs/pseudotime_mouse.executed.ipynb 
echo ALL_DONE
