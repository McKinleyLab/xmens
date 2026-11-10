#!/bin/bash
# Step 4 of 5: the tables and statistics behind the figures.
#   notebooks/4_statistics.ipynb (kernel rsc) writes the per-nucleus tables the figure notebooks read, the human x
#   mouse block of the joint graph, the share of cross-species edge weight per compartment, the mean total
#   cross-species edge weight per nucleus and SAMap's alignment scores between compartments, the agreement of the
#   transferred cell type with the human labels, and the permutation control of the pseudotime transfer (mouse
#   pseudotime shuffled across all mouse nuclei or within compartment, 1,000 times each, seed 20261022).
#   Reads:  $SAMAP_FINAL (step 3), $HUMAN_H5AD, $MOUSE_H5AD
#   Writes: $OUT/4_statistics/ tables/, graph_block/, edge_share/, alignment_scores.csv, label_agreement.csv, ptnull/,
#           statistics.json
#   Logs:   $OUT/notebooks/4_statistics.ipynb and logs/4_statistics_<job ID>.out
# Skipped when $OUT/4_statistics exists (the notebook writes 4_statistics.partial and renames it at the end); then it
# must have been made from $SAMAP_FINAL. The final objects are checked against the checksums config.sh pins.
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

# config.sh sets $OUT, $SAMAP_FINAL, $HUMAN_H5AD, $MOUSE_H5AD, run_notebook and exists.
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
    -p MOUSE_H5AD "$MOUSE_H5AD"
