#!/bin/bash
# Step 0b of 5: the phmmer hits -> the SAMap homology maps.
#   notebooks/0_homology.ipynb (kernel rsc) names each hit's proteins by their gene symbols and keeps, for each pair
#   of human and mouse genes, the best-scoring hit between their proteins (phmmer_to_samap_symbols.py of the run of
#   record).
#   Reads:  $PHMMER_HITS/{hu_vs_mu,mu_vs_hu}/<prefix>_<i>.tbl (step 0a), $HUMAN_PROTEOME, $MOUSE_PROTEOME
#   Writes: $OUT/0_homology/hsmm/hs_to_mm.txt and mm_to_hs.txt (the maps SAMap reads, in BLAST's tabular layout) and
#           homology.json (counts, known orthologs found, checksums)
#   Logs:   $OUT/notebooks/0_homology.ipynb and logs/0b_homology_<job ID>.out
# Skipped when both maps exist in $MAPS/hsmm (by default $OUT/0_homology/hsmm; then they must have been made from
# $PHMMER_HITS), or when a later step's input is supplied (SAMAP_JOINT or SAMAP_FINAL; see config.sh). The proteomes
# are checked against the checksums config.sh pins.
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

# config.sh sets $OUT, $PHMMER_HITS, $MAPS, the proteomes and their checksums, and the helper functions.
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
