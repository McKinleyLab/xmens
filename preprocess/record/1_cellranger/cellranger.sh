#!/bin/bash
#SBATCH -c 8
#SBATCH -t 5-00:00
#SBATCH -p eddy
#SBATCH --mem 100gb
#SBATCH --exclude=holygpu8a[26504-26506]
#SBATCH -o /n/eddy_lab/users/nhilgert/decidualization/revisions/velocity/%x/logs/%x_%j.out
#SBATCH -e /n/eddy_lab/users/nhilgert/decidualization/revisions/velocity/%x/logs/%x_%j.err

# Args passed from submit_cellranger.sh
LANE=$1
SAMPLE=$2
FASTQS=$3
REF=${4:-/n/eddy_lab/users/nhilgert/reference/transgenic_mus}

VELOCITY_DIR=/n/eddy_lab/users/nhilgert/decidualization/revisions/velocity
BASE=${VELOCITY_DIR}/${LANE}

echo "========================================="
echo "Cell Ranger count - ${LANE}"
echo "Sample: ${SAMPLE}"
echo "Fastqs: ${FASTQS}"
echo "Reference: ${REF}"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "========================================="

START_TIME=$SECONDS

cd ${BASE}

/n/eddy_lab/users/nhilgert/software/cellranger-7.2.0/cellranger count \
  --id=${LANE} \
  --transcriptome=${REF} \
  --fastqs=${FASTQS} \
  --sample=${SAMPLE} \
  --output-dir=${BASE}/cellranger \
  --localcores=8 \
  --localmem=96

EXIT_CODE=$?
ELAPSED=$(( SECONDS - START_TIME ))
HOURS=$(( ELAPSED / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
  echo "Cell Ranger count - ${LANE} COMPLETED SUCCESSFULLY"
else
  echo "Cell Ranger count - ${LANE} FAILED (exit code: $EXIT_CODE)"
  echo "Finished: $(date)"
  echo "Time taken: ${HOURS}h ${MINS}m ${SECS}s"
  exit 1
fi
echo "Finished: $(date)"
echo "Time taken: ${HOURS}h ${MINS}m ${SECS}s"
echo "========================================="
