#!/bin/bash
#SBATCH --job-name=phmmer_mu_vs_hu
#SBATCH --partition=eddy
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
#SBATCH --array=1-32
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/n/eddy_lab/users/nhilgert/data/proteomes/slurm_logs/phmmer_mu_vs_hu/array_%A_task_%a.out
#SBATCH --error=/n/eddy_lab/users/nhilgert/data/proteomes/slurm_logs/phmmer_mu_vs_hu/array_%A_task_%a.err

set -euo pipefail

PROT_DIR=/n/eddy_lab/users/nhilgert/data/proteomes
cd "${PROT_DIR}"
mkdir -p "${PROT_DIR}/slurm_logs/phmmer_mu_vs_hu"
mkdir -p "${PROT_DIR}/hm/mu_vs_hu"

CHUNK="${PROT_DIR}/chunks/mm/Mus_musculus.GRCm39.pep.all.fa.${SLURM_ARRAY_TASK_ID}"
TARGET="${PROT_DIR}/fasta/Homo_sapiens.GRCh38.pep.all.fa"
TBL="${PROT_DIR}/hm/mu_vs_hu/mu_${SLURM_ARRAY_TASK_ID}.tbl"
PHMMER_LOG="${PROT_DIR}/hm/mu_vs_hu/mu_${SLURM_ARRAY_TASK_ID}.log"

echo "=== TASK INFO ==="
echo "Array job:   ${SLURM_ARRAY_JOB_ID}"
echo "Task ID:     ${SLURM_ARRAY_TASK_ID} of 32"
echo "Host:        $(hostname)"
echo "Threads:     ${SLURM_CPUS_PER_TASK}"
echo "Query:       ${CHUNK} ($(grep -c '^>' ${CHUNK}) sequences)"
echo "Target:      ${TARGET} ($(grep -c '^>' ${TARGET}) sequences)"
echo "Output tbl:  ${TBL}"
echo "Start:       $(date)"
echo ""

time phmmer --cpu ${SLURM_CPUS_PER_TASK} -E 1e-5 \
    --tblout "${TBL}" \
    "${CHUNK}" "${TARGET}" > "${PHMMER_LOG}"

echo ""
echo "=== DONE ==="
echo "End:         $(date)"
echo "Hits in tbl: $(grep -cv '^#' ${TBL})"
