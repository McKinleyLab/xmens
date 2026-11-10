#!/bin/bash
# Compares the run in $OUT with the heatmaps of record (record/package/), for every step that has finished, so it can
# run at any point:
#   1 inputs   nuclei, pseudotime, gene lists and the counts of the drawn genes, against the package's inputs
#   2 monocle  Monocle's dispersion table and trend against the package's R port; fallback fits
#   3 figures  the z matrices, the MI gene order and every SVG against the figures of record
#   Writes: $OUT/compare_with_record_<date>_<time>.tsv, one row per comparison
#   Log:    logs/compare_with_record_<job ID>.out, which ends with the summary
#
# Usage:
#   sbatch compare_with_record.sbatch

#SBATCH --job-name=compare_with_record
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G  # one species' counts matrix (about 88 million non-zero counts)
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out  # job name and job ID; relative to the folder sbatch was run from

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# The pipeline folder holds config.sh. Run with bash, this script sits in it; under sbatch, SLURM runs a copy from its
# spool folder, and the pipeline folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the xmens_pt_heatmaps folder (no config.sh in $HERE)" >&2
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
    --record "$HERE/record" \
    --out "$OUT/compare_with_record_$(date +%Y%m%d_%H%M%S).tsv"
