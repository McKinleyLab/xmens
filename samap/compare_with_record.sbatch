#!/bin/bash
# Compares the run in $OUT with the SAMap run of record (record/), step by step, for every step that has finished,
# so it can run at any point:
#   0 homology    proteome chunks, phmmer hit tables and the homology maps against the run of record (24 Apr 2026)
#   1 inputs      counts, nuclei and genes against the inputs of record (inputs_raw/, 1 Sep 2026)
#   2 fit         alignment scores, the joint graph and UMAP against the fit of record (runs_raw/vanilla_raw)
#   3 transfer    Samap_Final.h5ad against the deposited one: obs column by column, every other dataset by checksum
#   4 statistics  tables, edge shares, label agreement and permutation statistics against the 3 Nov 2026 set
#   5 figures     against the 3 Nov 2026 set (record/4_5_figures/results_final/): every PNG pixel by pixel, every SVG
#                 apart from its drawing date and element ids, the CSV value by value
#   numbers       the numbers of record/numbers.csv that the run computes
#   Writes: $OUT/compare_with_record_<date>_<time>.csv, one row per comparison
#   Log:    logs/compare_with_record_<job ID>.out, which ends with the summary
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

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# The pipeline folder holds config.sh. Run with bash, this script sits in it; under sbatch, SLURM runs a copy from its
# spool folder, and the pipeline folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the samap folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, $PY and the start-from-earlier variables, whose files are compared as the run's.
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
