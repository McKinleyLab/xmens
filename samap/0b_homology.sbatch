#!/bin/bash
# Step 0b: the phmmer hits -> the homology table SAMap reads (notebooks/0_homology.ipynb).
# Reads $PHMMER_HITS (step 0a); writes $OUT/0_homology/hsmm/{hs_to_mm,mm_to_hs}.txt.
#
# Usage:
#   bash 0b_homology.sbatch --check
#   sbatch 0b_homology.sbatch

#SBATCH --job-name=0b_homology
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
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
    for direction in hu_vs_mu mu_vs_hu; do
        found=0
        if [ -d "$PHMMER_HITS/$direction" ]; then
            found=$(find "$PHMMER_HITS/$direction" -name '*.tbl' | wc -l)
        fi
        echo "reads:  $PHMMER_HITS/$direction/: $found of 32 hit tables"
    done
    echo "reads:  $HUMAN_PROTEOME ($(pinned "$HUMAN_PROTEOME" "$HUMAN_PROTEOME_SHA256"))"
    echo "        $MOUSE_PROTEOME ($(pinned "$MOUSE_PROTEOME" "$MOUSE_PROTEOME_SHA256"))"
    echo "writes: $MAPS/hsmm/hs_to_mm.txt ($(exists "$MAPS/hsmm/hs_to_mm.txt"))"
    echo "        $MAPS/hsmm/mm_to_hs.txt ($(exists "$MAPS/hsmm/mm_to_hs.txt"))"
    exit 0
fi

skip_if_unused 0b SAMAP_JOINT SAMAP_FINAL
if [ -e "$MAPS/hsmm/hs_to_mm.txt" ] && [ -e "$MAPS/hsmm/mm_to_hs.txt" ]; then
    if [ "$MAPS" = "$DEFAULT_MAPS" ]; then
        made_from "$MAPS/hsmm/homology.json" phmmer_hits "$PHMMER_HITS"
    fi
    echo "the homology maps exist in $MAPS/hsmm; skipped"
    exit 0
fi
if [ "$MAPS" != "$OUT/0_homology" ]; then
    echo "MAPS is set to $MAPS, which lacks the maps; step 0b writes only to $OUT/0_homology" >&2
    exit 1
fi

verify_sha256 "$HUMAN_PROTEOME" "$HUMAN_PROTEOME_SHA256"
verify_sha256 "$MOUSE_PROTEOME" "$MOUSE_PROTEOME_SHA256"
run_notebook "$HERE/notebooks/0_homology.ipynb" rsc "$OUT/notebooks/0_homology.ipynb" \
    -p OUT "$OUT" \
    -p PHMMER_HITS "$PHMMER_HITS" \
    -p HUMAN_PROTEOME "$HUMAN_PROTEOME" \
    -p MOUSE_PROTEOME "$MOUSE_PROTEOME"
