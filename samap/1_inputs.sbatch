#!/bin/bash
# Step 1: the final fibroblast/decidual objects -> the SAMap inputs (notebooks/1_inputs.ipynb).
# Reads $HUMAN_H5AD and $MOUSE_H5AD (config.sh); writes $OUT/1_inputs/.
#
# Usage:
#   bash 1_inputs.sbatch --check    print the inputs and outputs and whether they exist; run nothing
#   sbatch 1_inputs.sbatch

#SBATCH --job-name=1_inputs
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
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
    echo "reads:  $HUMAN_H5AD ($(pinned "$HUMAN_H5AD" "$HUMAN_H5AD_SHA256"))"
    echo "        $MOUSE_H5AD ($(pinned "$MOUSE_H5AD" "$MOUSE_H5AD_SHA256"))"
    echo "writes: $SAMAP_INPUTS/hs_fib_dec_all_raw.h5ad ($(exists "$SAMAP_INPUTS/hs_fib_dec_all_raw.h5ad"))"
    echo "        $SAMAP_INPUTS/mm_fib_dec_all_raw.h5ad ($(exists "$SAMAP_INPUTS/mm_fib_dec_all_raw.h5ad"))"
    exit 0
fi

skip_if_unused 1 SAMAP_JOINT SAMAP_FINAL
if [ -e "$SAMAP_INPUTS/hs_fib_dec_all_raw.h5ad" ] && [ -e "$SAMAP_INPUTS/mm_fib_dec_all_raw.h5ad" ]; then
    echo "the SAMap inputs exist in $SAMAP_INPUTS; skipped"
    exit 0
fi
if [ "$SAMAP_INPUTS" != "$OUT/1_inputs" ]; then
    echo "SAMAP_INPUTS is set to $SAMAP_INPUTS, which lacks the inputs; step 1 writes only to $OUT/1_inputs" >&2
    exit 1
fi

verify_sha256 "$HUMAN_H5AD" "$HUMAN_H5AD_SHA256"
verify_sha256 "$MOUSE_H5AD" "$MOUSE_H5AD_SHA256"
run_notebook "$HERE/notebooks/1_inputs.ipynb" samap "$OUT/notebooks/1_inputs.ipynb" \
    -p OUT "$OUT" \
    -p HUMAN_H5AD "$HUMAN_H5AD" \
    -p MOUSE_H5AD "$MOUSE_H5AD" \
    -p COUNTS_LAYER "$COUNTS_LAYER"
