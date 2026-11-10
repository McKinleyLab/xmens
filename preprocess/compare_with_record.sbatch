#!/bin/bash
# Compares the run in $OUT with the November 2026 run of record, sample by sample, for every step that has
# finished, so it can run at any point:
#   1 Cell Ranger  against the February 2026 run, byte for byte
#   2 CellBender   against the July 2026 run: cells called, and nuclei with identical counts
#   3 QC           against the record's QC logs (2026_11_06_qc_regen/logs)
#   4 and 5        the final objects and their positions against the regen objects
#   Writes: $OUT/compare_with_record_<date>_<time>.tsv, one row per sample
#   Log:    logs/compare_with_record_<job ID>.out, which ends with the summary
#
# Usage:
#   sbatch compare_with_record.sbatch

#SBATCH --job-name=compare_with_record
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G  # holds one species' record and run objects at once
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out  # job name and job ID; relative to the folder sbatch was run from

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

# config.sh sets $OUT and $PY.
source "$HERE/config.sh"

if [ ! -d "$OUT" ]; then
    echo "no run folder $OUT yet" >&2
    exit 1
fi

# The date and time in the table's name keep each comparison; none overwrites an earlier one.
"$PY" "$HERE/scripts/compare_with_record.py" \
    --run "$OUT" \
    --samples "$HERE/samples.tsv" \
    --out "$OUT/compare_with_record_$(date +%Y%m%d_%H%M%S).tsv"
