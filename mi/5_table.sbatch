#!/bin/bash
# Step 5 of 6: the supplementary table of pseudotime MI for one-to-one human-mouse orthologs, with a live overlap
# calculator on its front sheet (MI in bits).
#   Reads:  $OUT/inputs/<species>/pseudotime.npy (step 1, for the nuclei counts), $OUT/mi_mouse.csv and
#           $OUT/mi_human.csv (step 3), $OUT/overlap/ortholog_pairs_one2one.csv and methods_overlap.json (step 4)
#   Writes: $OUT/table/SuppTable_pseudotime_MI.xlsx, its data sheets as CSVs (SuppTable_pseudotime_MI__*.csv) and
#           SuppTable_pseudotime_MI_provenance.json
#   Logs:   $OUT/logs/5_table.log and logs/5_table_<job ID>.out
# Skipped if $OUT/table exists.
#
# Usage:
#   sbatch 5_table.sbatch

#SBATCH --job-name=5_table
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
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

if [ -e "$OUT/table" ]; then
    echo "$OUT/table exists; skipped"
    exit 0
fi
mkdir -p "$OUT/logs"
"$PY" "$HERE/scripts/supp_table.py" \
    --mi-mouse "$OUT/mi_mouse.csv" \
    --mi-human "$OUT/mi_human.csv" \
    --pairs "$OUT/overlap/ortholog_pairs_one2one.csv" \
    --overlap "$OUT/overlap/methods_overlap.json" \
    --inputs "$OUT/inputs" \
    --out-dir "$OUT/table" \
    2>&1 | tee "$OUT/logs/5_table.log"
