#!/bin/bash
# Usage: bash submit_cellranger.sh mck_11a [mck_12a ...]
# Submits a cellranger job for each lane, then chains a velocyto job on success.

VELOCITY_DIR=/n/eddy_lab/users/nhilgert/decidualization/revisions/velocity
LANES_DIR=/n/eddy_lab/users/nhilgert/decidualization/lanes
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for LANE in "$@"; do
    LANE_DIR=${VELOCITY_DIR}/${LANE}
    GEX_DIR=${LANES_DIR}/${LANE}/gex

    # Validate lane exists in lanes/
    if [ ! -d "${GEX_DIR}" ]; then
        echo "ERROR: No gex directory found for ${LANE} at ${GEX_DIR}, skipping."
        continue
    fi

    # Parse sample name from first fastq filename e.g. SI-TT-F1_S3_L004_R1_001.fastq.gz -> SI-TT-F1
    SAMPLE=$(ls "${GEX_DIR}"/*.fastq.gz 2>/dev/null | head -1 | xargs basename | grep -oP 'SI-[A-Z]+-[A-Z0-9]+')
    if [ -z "${SAMPLE}" ]; then
        echo "ERROR: Could not parse sample name from fastqs in ${GEX_DIR}, skipping."
        continue
    fi

    # Pre-create all output directories so SLURM log paths exist before job starts
    mkdir -p ${LANE_DIR}/{cellranger,velocyto,logs,scripts}

    # Submit cellranger, capture job ID
    CR_JOBID=$(sbatch --job-name=${LANE} \
        --parsable \
        ${SCRIPT_DIR}/cellranger.sh ${LANE} ${SAMPLE} ${GEX_DIR})
    echo "Submitted cellranger for ${LANE} (sample: ${SAMPLE}) — job ${CR_JOBID}"
    
    VE_JOBID=$(sbatch --job-name=${LANE} \
        --dependency=afterok:${CR_JOBID} \
        --parsable \
        ${SCRIPT_DIR}/velocyto.sh ${LANE})
    echo "Submitted velocyto  for ${LANE} — job ${VE_JOBID} (depends on ${CR_JOBID})"
done
