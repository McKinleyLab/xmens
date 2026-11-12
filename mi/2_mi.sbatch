#!/bin/bash
# Step 2 of 6: mutual information (MI) between each gene's expression and pseudotime, with 1,000 permutations.
# Each species' genes are split into chunks (config.sh: 200 mouse, 1,500 human) and each array task scores one chunk:
# tasks 0-199 are mouse chunks 0-199, tasks 200-1699 human chunks 0-1499.
#   scripts/mi_chunk.py: MI from scikit-learn's mutual_info_regression (3 neighbours, random_state 0) on all nuclei;
#     the chunk's 1,000 shuffles of the pseudotime are drawn once, seeded by the chunk number, and shared by its genes;
#     p = (shuffles with MI >= observed + 1) / 1,001.
#   Reads:  $OUT/inputs/<species>/ (step 1)
#   Writes: $OUT/chunks/<species>/mi_<chunk>.csv (gene, mi in nats, p)
#   Logs:   logs/2_mi_<job ID>_<task>.out
# In the run of record a mouse chunk (71-72 genes) took 3.7-6.1 h and 2.7 GB, a human chunk (10-11 genes) about 1 h.
# A task whose output exists is skipped, so after a failure the whole array can simply be submitted again.
#
# Usage:
#   bash 2_mi.sbatch --check                   count the finished chunks of each species; run nothing
#   sbatch 2_mi.sbatch                         all 1,700 tasks
#   sbatch --array=0,200 2_mi.sbatch           mouse chunk 0 and human chunk 0 only
#   SLURM_ARRAY_TASK_ID=200 bash 2_mi.sbatch   human chunk 0 in the current shell (an interactive session)

#SBATCH --job-name=2_mi
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00  # the run of record's limit for mouse chunks
#SBATCH --array=0-1699  # 200 mouse + 1,500 human chunks; override with sbatch --array
#SBATCH --output=logs/%x_%A_%a.out  # job name, job ID, task; relative to the folder sbatch was run from

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

# config.sh sets $OUT, $PY, $N_CHUNKS_MOUSE, $N_CHUNKS_HUMAN, n_chunks <species> and exists <path>.
source "$HERE/config.sh"

# --check: whether each species' inputs exist and how many of its chunks are done; then stop without running anything.
if [ "${1:-}" = "--check" ]; then
    for species in mouse human; do
        n_done=0
        if [ -d "$OUT/chunks/$species" ]; then
            n_done=$(find "$OUT/chunks/$species" -name 'mi_*.csv' | wc -l)
        fi
        echo "$species: inputs $(exists "$OUT/inputs/$species"); $n_done of $(n_chunks $species) chunks done"
    done
    exit 0
fi

# The array task ID, which sbatch sets, picks the species and chunk: mouse chunks first, then human.
TASK=${SLURM_ARRAY_TASK_ID:?submit with sbatch, or set SLURM_ARRAY_TASK_ID to a task number}
if [ "$TASK" -lt "$N_CHUNKS_MOUSE" ]; then
    SPECIES=mouse
    CHUNK=$TASK
elif [ "$TASK" -lt $((N_CHUNKS_MOUSE + N_CHUNKS_HUMAN)) ]; then
    SPECIES=human
    CHUNK=$((TASK - N_CHUNKS_MOUSE))
else
    echo "task $TASK is past the last chunk (tasks run 0 to $((N_CHUNKS_MOUSE + N_CHUNKS_HUMAN - 1)))" >&2
    exit 2
fi
RESULT=$OUT/chunks/$SPECIES/$(printf 'mi_%03d.csv' "$CHUNK")

if [ -e "$RESULT" ]; then
    echo "$RESULT exists; skipped"
    exit 0
fi
mkdir -p "$OUT/chunks/$SPECIES"
"$PY" "$HERE/scripts/mi_chunk.py" \
    --inputs "$OUT/inputs/$SPECIES" \
    --n-chunks "$(n_chunks $SPECIES)" \
    --chunk "$CHUNK" \
    --out "$RESULT"
