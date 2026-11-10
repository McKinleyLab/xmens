#!/bin/bash
# Step 3 of 6: per-chunk results -> one MI table per species, with Benjamini-Hochberg q over all of its genes.
#   Reads:  $CHUNKS (config.sh; by default step 2's $OUT/chunks/<species>/)
#   Writes: $OUT/mi_mouse.csv and $OUT/mi_human.csv (gene, mi in nats, p, q; sorted by q)
#   Logs:   $OUT/logs/3_combine_<species>.log and logs/3_combine_<job ID>.out
# Stops if a chunk is missing. A species whose table exists is skipped.
#
# Usage:
#   bash 3_combine.sbatch --check    print where the chunks are read from and how many there are; run nothing
#   sbatch 3_combine.sbatch

#SBATCH --job-name=3_combine
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
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

# config.sh sets $OUT, $PY, $CHUNKS, n_chunks <species> and exists <path>.
source "$HERE/config.sh"

if [ "${1:-}" = "--check" ]; then
    for species in mouse human; do
        chunks=${CHUNKS//\{species\}/$species}
        n_found=0
        if [ -d "$chunks" ]; then
            n_found=$(find "$chunks" -name 'mi_*.csv' | wc -l)
        fi
        echo "$species: $n_found of $(n_chunks $species) chunks in $chunks; table $(exists "$OUT/mi_$species.csv")"
    done
    exit 0
fi

mkdir -p "$OUT/logs"
for species in mouse human; do
    if [ -e "$OUT/mi_$species.csv" ]; then
        echo "$species: $OUT/mi_$species.csv exists; skipped"
        continue
    fi
    # The species' chunk folder: CHUNKS with {species} replaced (bash ${VAR//pattern/replacement}).
    "$PY" "$HERE/scripts/combine_chunks.py" \
        --chunks "${CHUNKS//\{species\}/$species}" \
        --n-chunks "$(n_chunks $species)" \
        --out "$OUT/mi_$species.csv" \
        2>&1 | tee "$OUT/logs/3_combine_$species.log"
done
