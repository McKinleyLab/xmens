#!/bin/bash
# Usage: bash submit_cellranger_human.sh mck_5a mck_5b mck_5c mck_6a mck_6b mck_6c
# Submits cellranger jobs against the human GRCh38 reference.

VELOCITY_DIR=/n/eddy_lab/users/nhilgert/decidualization/revisions/velocity
LANES_DIR=/n/eddy_lab/users/nhilgert/decidualization/lanes
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
HUMAN_REF=/n/eddy_lab/users/nhilgert/reference/human/refdata-gex-GRCh38-2024-A

for LANE in "$@"; do
    LANE_DIR=${VELOCITY_DIR}/${LANE}
    GEX_DIR=${LANES_DIR}/${LANE}/gex

    if [ ! -d "${GEX_DIR}" ]; then
        echo "ERROR: No gex directory found for ${LANE} at ${GEX_DIR}, skipping."
        continue
    fi

    SAMPLE=$(ls "${GEX_DIR}"/*.fastq.gz 2>/dev/null | head -1 | xargs basename | grep -oP 'SI-[A-Z]+-[A-Z0-9]+')
    if [ -z "${SAMPLE}" ]; then
        echo "ERROR: Could not parse sample name from fastqs in ${GEX_DIR}, skipping."
        continue
    fi

    # Wipe any previous (wrong) cellranger output before resubmitting
    rm -rf ${LANE_DIR}/cellranger
    mkdir -p ${LANE_DIR}/{cellranger,velocyto,logs,scripts}

    CR_JOBID=$(sbatch --job-name=${LANE} \
        --parsable \
        ${SCRIPT_DIR}/cellranger.sh ${LANE} ${SAMPLE} ${GEX_DIR} ${HUMAN_REF})
    echo "Submitted cellranger (HUMAN GRCh38) for ${LANE} (sample: ${SAMPLE}) — job ${CR_JOBID}"

    VE_JOBID=$(sbatch --job-name=${LANE} \
        --dependency=afterok:${CR_JOBID} \
        --parsable \
        ${SCRIPT_DIR}/velocyto.sh ${LANE})
    echo "Submitted velocyto for ${LANE} — job ${VE_JOBID} (depends on ${CR_JOBID})"
done
