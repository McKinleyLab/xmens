#!/bin/bash
# Compares the run in $OUT with the run of record (September-October 2026), for every step that has finished, so it
# can run at any point:
#   1 pooled objects   the loom's genes, nuclei and counts; the h5ad's nuclei, genes, raw counts and scVI latent
#   2 GRNBoost2        each run's sketch (its nuclei) and adjacencies (links and importances)
#   3 cisTarget        each run's regulons (names, targets, weights)
#   4 consensus        both floors (recurrence report, regulons, comparator AUCell) and the AUCell of record
#   5 correlations     rho per regulon and the Fig. 5E table, against the table of record
#   6 figures          Fig. 5E (the SVG, byte for byte) and Fig. S28 (its text and its two input tables)
#   Writes: $OUT/compare_with_record_<date>_<time>.tsv, one row per comparison
#   Log:    logs/compare_with_record_<job ID>.out, which ends with the summary
# The run-of-record files it reads are listed in record/SOURCES.tsv.
#
# Usage:
#   sbatch compare_with_record.sbatch

#SBATCH --job-name=compare_with_record
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G  # holds the record's and the run's pooled loom one block of nuclei at a time, and two AUCell tables
#SBATCH --time=04:00:00  # about 1 minute per GRNBoost2 run compared, for a full run of 200
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506  # GPU nodes
#SBATCH --output=logs/%x_%j.out  # job name and job ID; relative to the folder sbatch was run from

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# The pipeline folder holds config.sh and runs.tsv. Run with bash, this script sits in it; under sbatch, SLURM runs a
# copy from its spool folder, and the pipeline folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the scenic folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, $RECORD and $PY_SCENIC (the regulons are pickled pySCENIC objects, so the comparison runs in the
# pySCENIC environment).
source "$HERE/config.sh"

if [ ! -d "$OUT" ]; then
    echo "no run folder $OUT yet" >&2
    exit 1
fi

# The date and time in the table's name keep each comparison; none overwrites an earlier one.
"$PY_SCENIC" "$HERE/scripts/compare_with_record.py" \
    --run "$OUT" \
    --record "$RECORD" \
    --record-outputs "$HERE/record/outputs" \
    --runs "$HERE/runs.tsv" \
    --out "$OUT/compare_with_record_$(date +%Y%m%d_%H%M%S).tsv"
