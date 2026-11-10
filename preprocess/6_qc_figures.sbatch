#!/bin/bash
# Step 6: QC figures, one per sample, with the thresholds qc.py applied drawn on the distributions. Each array task
# runs one sample, task n being row n of samples.tsv. Changes no data; it reads step 3's objects and log.
#   Reads:  $OUT/raw/<sample>.h5ad, $OUT/qc/<sample>.h5ad and $OUT/logs/3_qc_<sample>.log (step 3)
#   Writes: $OUT/qc_figures/<sample>/<sample>_qc.png and summary.txt
#   Logs:   $OUT/logs/6_qc_figures_<sample>.log and logs/6_qc_figures_<job ID>_<task>.out
# The script recomputes the 3-MAD bounds as qc.py does and stops if they differ from the ones in the QC log.
# A sample whose figure exists is skipped.
#
# Usage:
#   sbatch 6_qc_figures.sbatch              all samples
#   sbatch --array=3 6_qc_figures.sbatch    row 3 only

#SBATCH --job-name=6_qc_figures
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
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

# config.sh sets $OUT, $PY, and field <row> <column>, which prints one value from samples.tsv.
source "$HERE/config.sh"

# The array task ID, which sbatch sets, is the row of samples.tsv this task runs; without one, the script stops
# with the message below.
ROW=${SLURM_ARRAY_TASK_ID:?submit with sbatch, or set SLURM_ARRAY_TASK_ID to a row of samples.tsv}
SAMPLE=$(field "$ROW" sample_id)
SPECIES=$(field "$ROW" species)
FIGURES=$OUT/qc_figures/$SAMPLE

if [ -e "$FIGURES/${SAMPLE}_qc.png" ]; then
    echo "$SAMPLE: $FIGURES/${SAMPLE}_qc.png exists; nothing to do"
    exit 0
fi
for input in "$OUT/raw/$SAMPLE.h5ad" "$OUT/qc/$SAMPLE.h5ad" "$OUT/logs/3_qc_$SAMPLE.log"; do
    if [ ! -e "$input" ]; then
        echo "$SAMPLE: input $input not found (run step 3 first)" >&2
        exit 1
    fi
done
mkdir -p "$FIGURES"

# Runs in rsc_uv ($PY), like ingest and concat; the figure needs no DoubletDetection.
"$PY" "$HERE/scripts/qc_report.py" \
    --raw "$OUT/raw/$SAMPLE.h5ad" \
    --qc "$OUT/qc/$SAMPLE.h5ad" \
    --qc-log "$OUT/logs/3_qc_$SAMPLE.log" \
    --sample "$SAMPLE" \
    --species "$SPECIES" \
    --outdir "$FIGURES" \
    --log "$OUT/logs/6_qc_figures_$SAMPLE.log"

# The log's last lines (bounds check, file written) also go to the SLURM log.
tail -n 3 "$OUT/logs/6_qc_figures_$SAMPLE.log"
