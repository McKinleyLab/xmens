#!/usr/bin/env bash
#SBATCH -p eddy                # partition
#SBATCH -c 1                   # number of CPU cores
#SBATCH --mem=16G              # memory per job
#SBATCH -t 10:00:00            # walltime
#SBATCH --array=0-999          # job array range
#SBATCH -o logs/mi_%A_%a.out   # standard output
#SBATCH -e logs/mi_%A_%a.err   # standard error
#SBATCH --job-name=MI          # job name
# (Optional) #SBATCH --exclude=holygpu8a2650[4-6]

set -euo pipefail

echo "Job ${SLURM_ARRAY_JOB_ID}/${SLURM_ARRAY_TASK_ID} started: $(date)"

# environment
module load python
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mi_env || { echo "conda env failed"; exit 1; }
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# paths
TOPDIR=/xmens/pseudotime/
H5AD=/xmens/pseudotime/data/20251008_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad

PY_SCRIPT=${TOPDIR}/scripts/mi_worker.py
OUTDIR=${TOPDIR}/results

cd "${TOPDIR}" || { echo "cd ${TOPDIR} failed"; exit 1; }

# ensure output dirs exist
mkdir -p logs "${OUTDIR}"

# calculate n chunks from the array range (1000 in this case)
NCHUNKS=$(( SLURM_ARRAY_TASK_MAX - SLURM_ARRAY_TASK_MIN + 1 ))

# ---------- launch ----------
python "${PY_SCRIPT}" \
  --h5ad "${H5AD}" \
  --chunk "${SLURM_ARRAY_TASK_ID}" \
  --n-chunks "${NCHUNKS}" \
  --outdir "${OUTDIR}" \
  --perms 1001

echo "Task ${SLURM_ARRAY_TASK_ID} finished: $(date)"