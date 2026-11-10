#!/bin/bash
#SBATCH --partition=eddy
#SBATCH --time=0-14:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --account=eddy_lab
#SBATCH --array=0-199%100
#SBATCH --job-name=mi_pt_harmonized_mouse_nall
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
#SBATCH --output=/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/logs/mi_dpt/pt_harmonized_mouse_nall_%a.out
export XDG_CACHE_HOME=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch NUMBA_CACHE_DIR=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch/numba MPLCONFIGDIR=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch/mpl TMPDIR=/n/lab_storage/eddy_lab/users/nhilgert/tmp
/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/workflow/scripts/mi/mi_pseudotime.py worker --outdir /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/results/mi_dpt/pt_harmonized_mouse_nall --chunk $SLURM_ARRAY_TASK_ID --perms 1001 --subsample-cells 0
