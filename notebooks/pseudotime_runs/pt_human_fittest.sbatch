#!/bin/bash
#SBATCH -J pt_human_fittest
#SBATCH -p eddy
#SBATCH --account=eddy_lab
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 03:00:00
#SBATCH -o /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks/pseudotime_runs/pt_human_fittest-%j.out
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
set -euo pipefail
C=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch
export HOME=/n/lab_storage/eddy_lab/users/nhilgert/tmp MPLCONFIGDIR=$C/mpl NUMBA_CACHE_DIR=$C/numba IPYTHONDIR=$C/ipython XDG_CACHE_HOME=$C \
       JUPYTER_DATA_DIR=$C/jupyter_data JUPYTER_CONFIG_DIR=$C/jupyter_config JUPYTER_RUNTIME_DIR=/n/lab_storage/eddy_lab/users/nhilgert/tmp/jupyter_runtime_pt JAX_PLATFORMS=cpu
mkdir -p $JUPYTER_RUNTIME_DIR
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks
hostname; lscpu | grep 'Model name'; nproc
/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python -m papermill --log-output --kernel python3 pseudotime_human.ipynb /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks/pseudotime_runs/pseudotime_human.fit_latent_test_2_epochs.executed.ipynb -p FIT_LATENT True -p FIT_MAX_EPOCHS 2 -p OUTPUT_OBJECT /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_pseudotime_notebooks/test_fit_latent_2_epochs/Cevrim_Human_Only_Fib_Dec_Cells_pseudotime_TEST.h5ad
echo ALL_DONE
