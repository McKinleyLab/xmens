#!/bin/bash
# Step 4: the tables and statistics behind the figures (notebooks/4_statistics.ipynb).
# Reads $SAMAP_FINAL (step 3) and the final objects; writes $OUT/4_statistics/.
#
# Usage:
#   bash 4_statistics.sbatch --check
#   sbatch 4_statistics.sbatch

#SBATCH --job-name=4_statistics
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
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

if [ "${1:-}" = "--check" ]; then
    echo "reads:  $SAMAP_FINAL ($(exists "$SAMAP_FINAL"))"
    echo "        $HUMAN_H5AD ($(pinned "$HUMAN_H5AD" "$HUMAN_H5AD_SHA256"))"
    echo "        $MOUSE_H5AD ($(pinned "$MOUSE_H5AD" "$MOUSE_H5AD_SHA256"))"
    echo "writes: $OUT/4_statistics ($(exists "$OUT/4_statistics"))"
    exit 0
fi

if [ -e "$OUT/4_statistics" ]; then
    made_from "$OUT/4_statistics/statistics.json" samap_final "$SAMAP_FINAL"
    echo "$OUT/4_statistics exists; skipped"
    exit 0
fi

verify_sha256 "$HUMAN_H5AD" "$HUMAN_H5AD_SHA256"
verify_sha256 "$MOUSE_H5AD" "$MOUSE_H5AD_SHA256"
run_notebook "$HERE/notebooks/4_statistics.ipynb" rsc "$OUT/notebooks/4_statistics.ipynb" \
    -p OUT "$OUT" \
    -p SAMAP_FINAL "$SAMAP_FINAL" \
    -p HUMAN_H5AD "$HUMAN_H5AD" \
    -p MOUSE_H5AD "$MOUSE_H5AD" \
    -p X_COLUMN "$X_COLUMN" -p Y_COLUMN "$Y_COLUMN"
