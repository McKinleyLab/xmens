#!/bin/bash
# Step 0a: phmmer (HMMER 3.4) searches of each proteome against the other, in 64 array tasks: tasks 1-32 are human
# chunks against mouse, tasks 33-64 mouse chunks against human. Writes $OUT/0_homology/phmmer/.
# The tasks that hold titin (12-23, 42, 56, 57) fail at 8G: sbatch --array=12-23,42,56,57 --mem=64G 0a_phmmer.sbatch
#
# Usage:
#   bash 0a_phmmer.sbatch --check             how many of the 64 hit tables exist; run nothing
#   sbatch 0a_phmmer.sbatch                   all 64 tasks
#   sbatch --array=1,33 0a_phmmer.sbatch      the first chunk of each direction only

#SBATCH --job-name=0a_phmmer
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=4  # phmmer --cpu
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --array=1-64  # 32 chunks per direction
#SBATCH --output=logs/%x_%A_%a.out  # job name, job ID, task; relative to the folder sbatch was run from

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

N_CHUNKS=32
HOMOLOGY=$OUT/0_homology

# Paths for array task $1.
task_paths() {
    if [ "$1" -le "$N_CHUNKS" ]; then
        DIRECTION=hu_vs_mu PREFIX=hu SPECIES=hs QUERY=$HUMAN_PROTEOME TARGET=$MOUSE_PROTEOME I=$1
    else
        DIRECTION=mu_vs_hu PREFIX=mu SPECIES=mm QUERY=$MOUSE_PROTEOME TARGET=$HUMAN_PROTEOME I=$(($1 - N_CHUNKS))
    fi
    CHUNK=$HOMOLOGY/chunks/$SPECIES/$(basename "$QUERY").$I
    TBL=$HOMOLOGY/phmmer/$DIRECTION/${PREFIX}_$I.tbl
    PHMMER_LOG=$HOMOLOGY/phmmer/$DIRECTION/${PREFIX}_$I.log
}

if [ "${1:-}" = "--check" ]; then
    echo "reads:  $HUMAN_PROTEOME ($(pinned "$HUMAN_PROTEOME" "$HUMAN_PROTEOME_SHA256"))"
    echo "        $MOUSE_PROTEOME ($(pinned "$MOUSE_PROTEOME" "$MOUSE_PROTEOME_SHA256"))"
    echo "        $PHMMER ($(exists "$PHMMER"))"
    for direction in hu_vs_mu mu_vs_hu; do
        found=0
        if [ -d "$HOMOLOGY/phmmer/$direction" ]; then
            found=$(find "$HOMOLOGY/phmmer/$direction" -name '*.tbl' | wc -l)
        fi
        echo "writes: $HOMOLOGY/phmmer/$direction/: $found of $N_CHUNKS hit tables"
    done
    exit 0
fi

skip_if_unused 0a PHMMER_HITS MAPS SAMAP_JOINT SAMAP_FINAL

TASK=${SLURM_ARRAY_TASK_ID:?submit with sbatch, or set SLURM_ARRAY_TASK_ID to a task, 1-64}
task_paths "$TASK"
if [ -e "$TBL" ]; then
    echo "$TBL exists; skipped"
    exit 0
fi
verify_sha256 "$HUMAN_PROTEOME" "$HUMAN_PROTEOME_SHA256"
verify_sha256 "$MOUSE_PROTEOME" "$MOUSE_PROTEOME_SHA256"
mkdir -p "$(dirname "$CHUNK")" "$(dirname "$TBL")"

# The chunk: proteins I, I+32, I+64, ... of the query proteome.
awk -v n="$N_CHUNKS" -v i="$I" '/^>/ { k++ } k > 0 && (k - 1) % n == i - 1' "$QUERY" > "$CHUNK.partial"
mv "$CHUNK.partial" "$CHUNK"

# --noali leaves the alignments out of the log; the hit table is the same either way.
"$PHMMER" --cpu "${SLURM_CPUS_PER_TASK:-4}" -E 1e-5 --noali --tblout "$TBL.partial" "$CHUNK" "$TARGET" > "$PHMMER_LOG"
mv "$TBL.partial" "$TBL"
echo "$TBL: $(grep -cv '^#' "$TBL") hits for $(grep -c '^>' "$CHUNK") query proteins"
