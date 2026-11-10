#!/bin/bash
# Compares the run in $OUT with the run of record, for every step that has finished, so it can run at any point:
#   1 inputs   against the record's expr.npz, genes.npy and pseudotime.npy (xmens_final/results/mi_dpt/)
#   2 chunks   every finished chunk against the record's results/mi_<chunk>.csv
#   3 combine  against the record's mi_all.csv
#   4 overlap  against xmens_final/results/mi_dpt/overlap_of_record/
#   5 table    against the approved workbook, record/5_table/SuppTable_pseudotime_MI.xlsx
#   6 fig5f    against overlap_of_record/crossspecies_heatmap_selection/ and the 28 genes drawn (record/6_fig5f/)
#   Writes: $OUT/compare_with_record_<date>_<time>.tsv, one row per comparison
#   Log:    logs/compare_with_record_<job ID>.out, which ends with the summary
#
# Usage:
#   sbatch compare_with_record.sbatch

#SBATCH --job-name=compare_with_record
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G  # one species' expression matrix from the run and from the record at once (6 GB)
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
    echo "submit this from the mi folder (no config.sh in $HERE)" >&2
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
    --record-table "$HERE/record/5_table/SuppTable_pseudotime_MI.xlsx" \
    --panel-genes "$HERE/record/6_fig5f/fig5f_panel_genes.csv" \
    --out "$OUT/compare_with_record_$(date +%Y%m%d_%H%M%S).tsv"
