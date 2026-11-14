#!/bin/bash
# Step 0c: renames genes in the homology table to the names they have in the objects, matched by Ensembl gene ID
# (scripts/relabel_maps.py). Reads $SOURCE_MAPS/hsmm/; writes $OUT/0_homology_relabeled/.
# Then run steps 2-5 with MAPS=$OUT/0_homology_relabeled.
#
# Usage:
#   export OUT=... SAMAP_INPUTS=... SOURCE_MAPS=$OUT/0_homology; sbatch 0c_relabel.sbatch

#SBATCH --job-name=0c_relabel
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
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
SOURCE_MAPS=${SOURCE_MAPS:-$NH/decidualization/revisions/samap_integrate/homology/hsmm_hmmer}
DEST=$OUT/0_homology_relabeled

if [ -e "$DEST" ]; then
    echo "$DEST exists; skipped"
    exit 0
fi
mkdir -p "$OUT/logs"
"$PY" "$HERE/scripts/relabel_maps.py" \
    --maps "$SOURCE_MAPS" \
    --human-proteome "$HUMAN_PROTEOME" --human-proteome-sha256 "$HUMAN_PROTEOME_SHA256" \
    --mouse-proteome "$MOUSE_PROTEOME" --mouse-proteome-sha256 "$MOUSE_PROTEOME_SHA256" \
    --human-h5ad "$SAMAP_INPUTS/hs_fib_dec_all_raw.h5ad" \
    --mouse-h5ad "$SAMAP_INPUTS/mm_fib_dec_all_raw.h5ad" \
    --out "$DEST" \
    2>&1 | tee "$OUT/logs/0c_relabel.log"
