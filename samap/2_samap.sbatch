#!/bin/bash
# Step 2: the SAMap fit (notebooks/2_samap.ipynb), about 1 h.
# Reads $SAMAP_INPUTS (step 1) and $MAPS/hsmm/ (step 0); writes $OUT/2_samap/.
#
# Usage:
#   bash 2_samap.sbatch --check
#   sbatch 2_samap.sbatch

#SBATCH --job-name=2_samap
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=8  # SAMap's ncpus
#SBATCH --mem=128G
#SBATCH --time=06:00:00
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
    for f in "$SAMAP_INPUTS/hs_fib_dec_all_raw.h5ad" "$SAMAP_INPUTS/mm_fib_dec_all_raw.h5ad" \
             "$MAPS/hsmm/hs_to_mm.txt" "$MAPS/hsmm/mm_to_hs.txt"; do
        echo "reads:  $f ($(exists "$f"))"
    done
    echo "writes: $SAMAP_JOINT ($(exists "$SAMAP_JOINT")) and the rest of $OUT/2_samap/"
    exit 0
fi

skip_if_unused 2 SAMAP_FINAL
if [ -e "$SAMAP_JOINT" ]; then
    if [ "$SAMAP_JOINT" = "$DEFAULT_SAMAP_JOINT" ]; then
        made_from "$OUT/2_samap/meta.json" inputs/hs/path "$SAMAP_INPUTS/hs_fib_dec_all_raw.h5ad"
        made_from "$OUT/2_samap/meta.json" maps/hs_to_mm.txt/path "$MAPS/hsmm/hs_to_mm.txt"
    fi
    echo "$SAMAP_JOINT exists; skipped"
    exit 0
fi
if [ "$SAMAP_JOINT" != "$OUT/2_samap/samap_joint.h5ad" ]; then
    echo "SAMAP_JOINT is set to $SAMAP_JOINT, which does not exist; step 2 writes only to $OUT/2_samap" >&2
    exit 1
fi

run_notebook "$HERE/notebooks/2_samap.ipynb" samap "$OUT/notebooks/2_samap.ipynb" \
    -p OUT "$OUT" \
    -p INPUTS "$SAMAP_INPUTS" \
    -p MAPS "$MAPS" \
    -p NCPUS "${SLURM_CPUS_PER_TASK:-8}"
