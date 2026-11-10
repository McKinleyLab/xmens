#!/bin/bash
# Step 1 of 5: the final fibroblast/decidual objects -> the SAMap inputs.
#   notebooks/1_inputs.ipynb (kernel samap) puts each object's raw counts (layer counts) in X, keeps the labels in
#   obs and drops the layers, embeddings and .raw; SAM normalizes the counts itself in step 2.
#   Reads:  $HUMAN_H5AD, $MOUSE_H5AD (config.sh)
#   Writes: $OUT/1_inputs/hs_fib_dec_all_raw.h5ad, mm_fib_dec_all_raw.h5ad and inputs.json (counts and checksums)
#   Logs:   $OUT/notebooks/1_inputs.ipynb (the executed notebook) and logs/1_inputs_<job ID>.out
# Skipped when both inputs exist in $SAMAP_INPUTS (by default $OUT/1_inputs), or when a later step's input is
# supplied (SAMAP_JOINT or SAMAP_FINAL; see config.sh). The final objects are checked against the checksums config.sh
# pins.
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

# config.sh sets $OUT, $SAMAP_INPUTS, the final objects and their checksums, and the helper functions.
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
    -p MOUSE_H5AD "$MOUSE_H5AD"
