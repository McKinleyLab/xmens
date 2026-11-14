#!/bin/bash
# Compares the run in $OUT with the SAMap run of record (record/), step by step, for every step that has finished.
# Writes $OUT/compare_with_record_<date>_<time>.csv.
#
# Usage:
#   sbatch compare_with_record.sbatch

#SBATCH --job-name=compare_with_record
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G  # holds the record's and the run's joint objects at once
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out  # job name and job ID; relative to the folder sbatch was run from

set -euo pipefail

# The pipeline folder: this script's folder, or under sbatch the folder the job was submitted from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the samap folder (no config.sh in $HERE)" >&2
    exit 2
fi

source "$HERE/config.sh"

if [ ! -d "$OUT" ]; then
    echo "no run folder $OUT yet" >&2
    exit 1
fi

# The date and time in the table's name keep each comparison; none overwrites an earlier one.
"$PY" "$HERE/scripts/compare_with_record.py" \
    --run "$OUT" \
    --phmmer-hits "$PHMMER_HITS" \
    --maps "$MAPS" \
    --inputs "$SAMAP_INPUTS" \
    --joint "$SAMAP_JOINT" \
    --samap-final "$SAMAP_FINAL" \
    --record "$HERE/record" \
    --out "$OUT/compare_with_record_$(date +%Y%m%d_%H%M%S).csv"
