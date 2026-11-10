#!/bin/bash
# Step 3 of 5: mouse cell type and pseudotime transferred onto the human nuclei, and Samap_Final.h5ad.
#   notebooks/3_transfer.ipynb (kernel rsc) transfers both over the human x mouse block of the joint SAMap graph and
#   writes the deposited object: the joint object with obs cut to species, sample, timepoint, cell_type,
#   transferred_mouse_pseudotime and transferred_mouse_fibdec.
#   Reads:  $SAMAP_JOINT (step 2), $HUMAN_H5AD, $MOUSE_H5AD (labels and the mouse pseudotime of record)
#   Writes: $OUT/3_transfer/ Samap_Final.h5ad, transfers.csv (one row per human nucleus), summary.json
#   Logs:   $OUT/notebooks/3_transfer.ipynb and logs/3_transfer_<job ID>.out
# Skipped when $SAMAP_FINAL exists (by default $OUT/3_transfer/Samap_Final.h5ad; then it must have been made from
# $SAMAP_JOINT). The final objects are checked against the checksums config.sh pins.
#
# Usage:
#   bash 3_transfer.sbatch --check
#   sbatch 3_transfer.sbatch

#SBATCH --job-name=3_transfer
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
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
    echo "submit this from the samap folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, $SAMAP_JOINT, $SAMAP_FINAL, $HUMAN_H5AD, $MOUSE_H5AD, run_notebook and exists.
source "$HERE/config.sh"

if [ "${1:-}" = "--check" ]; then
    echo "reads:  $SAMAP_JOINT ($(exists "$SAMAP_JOINT"))"
    echo "        $HUMAN_H5AD ($(pinned "$HUMAN_H5AD" "$HUMAN_H5AD_SHA256"))"
    echo "        $MOUSE_H5AD ($(pinned "$MOUSE_H5AD" "$MOUSE_H5AD_SHA256"))"
    echo "writes: $SAMAP_FINAL ($(exists "$SAMAP_FINAL")) and the rest of $OUT/3_transfer/"
    exit 0
fi

if [ -e "$SAMAP_FINAL" ]; then
    if [ "$SAMAP_FINAL" = "$DEFAULT_SAMAP_FINAL" ]; then
        made_from "$OUT/3_transfer/summary.json" samap_final/joint "$SAMAP_JOINT"
    fi
    echo "$SAMAP_FINAL exists; skipped"
    exit 0
fi
if [ "$SAMAP_FINAL" != "$OUT/3_transfer/Samap_Final.h5ad" ]; then
    echo "SAMAP_FINAL is set to $SAMAP_FINAL, which does not exist; step 3 writes only to $OUT/3_transfer" >&2
    exit 1
fi

verify_sha256 "$HUMAN_H5AD" "$HUMAN_H5AD_SHA256"
verify_sha256 "$MOUSE_H5AD" "$MOUSE_H5AD_SHA256"
run_notebook "$HERE/notebooks/3_transfer.ipynb" rsc "$OUT/notebooks/3_transfer.ipynb" \
    -p OUT "$OUT" \
    -p JOINT "$SAMAP_JOINT" \
    -p HUMAN_H5AD "$HUMAN_H5AD" \
    -p MOUSE_H5AD "$MOUSE_H5AD"
