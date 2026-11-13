#!/bin/bash
# Step 4 of 6: overlap of the human and mouse top MI genes among one-to-one orthologs: the Methods numbers (N 9,587,
# K 1,133, n 1,201, k 373, p = 1.95 x 10^-84 in the run of record) and the percentile sweep of Supplementary Fig. S29C-D.
#   Reads:  $OUT/mi_mouse.csv and $OUT/mi_human.csv (step 3), $ORTHOLOGS, and the gene IDs (var) of $H5AD_MOUSE and
#           $H5AD_HUMAN (config.sh), to match the pairs to the MI tables by Ensembl gene ID
#   Writes: $OUT/overlap/ortholog_pairs_one2one.csv, methods_overlap.json and overlap_percentile_sweep.csv
#   Logs:   $OUT/logs/4_overlap.log and logs/4_overlap_<job ID>.out
# Skipped if $OUT/overlap exists.
#
# Usage:
#   sbatch 4_overlap.sbatch

#SBATCH --job-name=4_overlap
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
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

# config.sh sets $OUT, $PY, $ORTHOLOGS, $H5AD_MOUSE and $H5AD_HUMAN.
source "$HERE/config.sh"

if [ -e "$OUT/overlap" ]; then
    echo "$OUT/overlap exists; skipped"
    exit 0
fi
mkdir -p "$OUT/logs"
"$PY" "$HERE/scripts/ortholog_overlap.py" \
    --mi-mouse "$OUT/mi_mouse.csv" \
    --mi-human "$OUT/mi_human.csv" \
    --orthologs "$ORTHOLOGS" \
    --h5ad-mouse "$H5AD_MOUSE" \
    --h5ad-human "$H5AD_HUMAN" \
    --out-dir "$OUT/overlap" \
    2>&1 | tee "$OUT/logs/4_overlap.log"
