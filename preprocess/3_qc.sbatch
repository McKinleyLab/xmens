#!/bin/bash
# Step 3 of 5: CellBender counts -> per-sample QC'd AnnData. Each array task runs one sample, task n being row n of
# samples.tsv.
#   scripts/ingest.py reads $CELLBENDER_H5 (config.sh; by default step 2's <sample>_cellbender_filtered.h5),
#     attaches Nick's Slide-tags positions ($COORDS_DIR/<lane>_coords.csv) and writes $OUT/raw/<sample>.h5ad with
#     every CellBender cell.
#   scripts/qc.py applies the thresholds (at least 500 UMIs; 3 MADs on log10 total counts and on log10 genes; at
#     most 8% mitochondrial) and DoubletDetection 4.2, and writes $OUT/qc/<sample>.h5ad.
#   Logs: $OUT/logs/3_ingest_<sample>.log and 3_qc_<sample>.log (every QC number), logs/3_qc_<job ID>_<task>.out
# Resources are those of the November 2026 run of record. A step whose output exists is skipped.
#
# Usage:
#   bash 3_qc.sbatch --check                  print each sample's inputs and outputs; run nothing
#   sbatch 3_qc.sbatch                        all samples
#   sbatch --array=3 3_qc.sbatch              row 3 only
#   SLURM_ARRAY_TASK_ID=3 bash 3_qc.sbatch    row 3 in the current shell (an interactive session)

#SBATCH --job-name=3_qc
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --constraint="genoa|turin"  # AMD nodes for reproducibility
#SBATCH --cpus-per-task=8  # qc.py runs DoubletDetection on all cores (n_jobs=-1)
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=1-15  # one task per row of samples.tsv
#SBATCH --output=logs/%x_%A_%a.out  # job name, job ID, task; relative to the folder sbatch was run from

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# The pipeline folder holds config.sh and samples.tsv. Run with bash, this script sits in it; under sbatch, SLURM
# runs a copy from its spool folder, and the pipeline folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the preprocess folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, $CELLBENDER_H5, $COORDS_DIR, $PY, $PY_QC, N_SAMPLES, and field <row> <column>, which prints
# one value from samples.tsv.
source "$HERE/config.sh"

# Sets SAMPLE, LANE, SPECIES, COUNTS, COORDS, RAW_H5AD and QC_H5AD for row $1 of samples.tsv. COUNTS is
# CELLBENDER_H5 with {sample} replaced by the row's sample (bash ${VAR//pattern/replacement}).
sample_paths() {
    SAMPLE=$(field "$1" sample_id)
    LANE=$(field "$1" lane)
    SPECIES=$(field "$1" species)
    COUNTS=${CELLBENDER_H5//\{sample\}/$SAMPLE}
    COORDS=$COORDS_DIR/${LANE}_coords.csv
    RAW_H5AD=$OUT/raw/$SAMPLE.h5ad
    QC_H5AD=$OUT/qc/$SAMPLE.h5ad
}

# Prints "found" or "missing", for --check.
exists() {
    if [ -e "$1" ]; then
        echo found
    else
        echo missing
    fi
}

# --check: for every row, print whether the inputs and outputs exist and where the counts come from; then stop
# without running anything.
if [ "${1:-}" = "--check" ]; then
    for row in $(seq 1 "$N_SAMPLES"); do
        sample_paths "$row"
        echo "$row $SAMPLE ($SPECIES, $LANE): counts $(exists "$COUNTS"), coords $(exists "$COORDS")," \
            "raw $(exists "$RAW_H5AD"), qc $(exists "$QC_H5AD")"
        echo "    counts: $COUNTS"
    done
    exit 0
fi

# The array task ID, which sbatch sets, is the row of samples.tsv this task runs; without one, the script stops
# with the message below.
ROW=${SLURM_ARRAY_TASK_ID:?submit with sbatch, or set SLURM_ARRAY_TASK_ID to a row of samples.tsv}
sample_paths "$ROW"
mkdir -p "$OUT/raw" "$OUT/qc" "$OUT/logs"

# Ingest runs in rsc_uv ($PY). Its output is kept, so a resubmitted task goes straight to QC.
if [ ! -e "$RAW_H5AD" ]; then
    "$PY" "$HERE/scripts/ingest.py" \
        --counts "$COUNTS" \
        --coords "$COORDS" \
        --sample "$SAMPLE" \
        --species "$SPECIES" \
        --timepoint "$(field "$ROW" timepoint)" \
        --lane "$LANE" \
        --modality "$(field "$ROW" modality)" \
        --out "$RAW_H5AD" \
        --log "$OUT/logs/3_ingest_$SAMPLE.log"
fi

# QC runs in the DoubletDetection 4.2 environment ($PY_QC); its thresholds and settings are qc.py's defaults.
if [ ! -e "$QC_H5AD" ]; then
    "$PY_QC" "$HERE/scripts/qc.py" \
        --in "$RAW_H5AD" \
        --out "$QC_H5AD" \
        --sample "$SAMPLE" \
        --species "$SPECIES" \
        --log "$OUT/logs/3_qc_$SAMPLE.log"
fi

# The QC log's last lines (doublets, nuclei kept, file written) also go to the SLURM log.
tail -n 3 "$OUT/logs/3_qc_$SAMPLE.log"
