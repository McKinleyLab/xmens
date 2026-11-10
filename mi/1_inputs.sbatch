#!/bin/bash
# Step 1 of 6: final objects -> the expression matrix, pseudotime and gene list that step 2 scores, per species.
#   scripts/mi_inputs.py keeps the nuclei with a pseudotime and the genes detected in at least 1% of them with at least
#     10 counts in total (14,261 mouse and 15,008 human genes in the run of record), and computes log1p(CP10K) from the
#     raw counts of those genes.
#   Reads:  $H5AD_MOUSE and $H5AD_HUMAN (config.sh)
#   Writes: $OUT/inputs/<species>/expression.npz, pseudotime.npy and genes.txt
#   Logs:   $OUT/logs/1_inputs_<species>.log and logs/1_inputs_<job ID>.out
# A species whose inputs folder exists is skipped; the folder is written under a temporary name and renamed when
# complete.
#
# Usage:
#   bash 1_inputs.sbatch --check    print the inputs and outputs and whether they exist; run nothing
#   sbatch 1_inputs.sbatch

#SBATCH --job-name=1_inputs
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G  # an object and a dense float64 copy of its expression; the test run took 2 min and 5.5 GB
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

# config.sh sets $OUT, $PY, $H5AD_MOUSE, $H5AD_HUMAN, h5ad <species> and exists <path>.
source "$HERE/config.sh"

if [ "${1:-}" = "--check" ]; then
    for species in mouse human; do
        echo "$species: object $(exists "$(h5ad $species)"), inputs $(exists "$OUT/inputs/$species")"
        echo "    object: $(h5ad $species)"
    done
    exit 0
fi

mkdir -p "$OUT/logs"
for species in mouse human; do
    if [ -e "$OUT/inputs/$species" ]; then
        echo "$species: $OUT/inputs/$species exists; skipped"
        continue
    fi
    "$PY" "$HERE/scripts/mi_inputs.py" \
        --h5ad "$(h5ad $species)" \
        --out-dir "$OUT/inputs/$species" \
        2>&1 | tee "$OUT/logs/1_inputs_$species.log"
done
