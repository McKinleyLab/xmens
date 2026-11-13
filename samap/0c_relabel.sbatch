#!/bin/bash
# Step 0c: the homology maps with each gene under the name it has in the objects (scripts/relabel_maps.py). Step 0b
# names genes by the Ensembl 113 proteomes' gene_symbol; genes renamed since the objects' annotation (mouse Sept2, now
# Septin2) are renamed back to the objects' names, matched by Ensembl gene ID, so that SAMap finds their homology.
#   Reads:  $SOURCE_MAPS/hsmm/{hs_to_mm,mm_to_hs}.txt (default: the maps of record), $HUMAN_PROTEOME and
#           $MOUSE_PROTEOME (config.sh, checked against their sha256), $SAMAP_INPUTS/{hs,mm}_fib_dec_all_raw.h5ad
#   Writes: $OUT/0_homology_relabeled/ hsmm/hs_to_mm.txt, hsmm/mm_to_hs.txt, renamed_genes.tsv, not_renamed.tsv,
#           relabel.json
#   Logs:   $OUT/logs/0c_relabel.log and logs/0c_relabel_<job ID>.out
# Skipped if $OUT/0_homology_relabeled exists. Then run steps 2-5 with MAPS=$OUT/0_homology_relabeled.
#
# Usage:
#   export OUT=... SAMAP_INPUTS=...; sbatch 0c_relabel.sbatch

#SBATCH --job-name=0c_relabel
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
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
    echo "submit this from the samap folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, $PY, $SAMAP_INPUTS and the proteomes with their sha256.
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
