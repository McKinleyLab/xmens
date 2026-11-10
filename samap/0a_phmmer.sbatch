#!/bin/bash
# Step 0a of 5: phmmer (HMMER 3.4) searches of each proteome against the other, the raw material of the SAMap
# homology maps. Each array task searches one chunk of query proteins against the whole other proteome: tasks 1-32
# human chunks against the mouse proteins (hu_vs_mu), tasks 33-64 mouse chunks against the human proteins (mu_vs_hu).
# Chunk i of a proteome holds its proteins i, i+32, i+64, ... (round robin), which reproduces the chunks of the run
# of record byte for byte.
#   Reads:  $HUMAN_PROTEOME, $MOUSE_PROTEOME, $PHMMER (config.sh)
#   Writes: $OUT/0_homology/chunks/{hs,mm}/<proteome>.<i>, and $OUT/0_homology/phmmer/hu_vs_mu/hu_<i>.tbl or
#           mu_vs_hu/mu_<i>.tbl (phmmer --tblout: one line per hit, read by step 0b) with phmmer's summary output
#           beside it (<prefix>_<i>.log)
#   Logs:   logs/0a_phmmer_<job ID>_<task>.out
# Resources are those of the run of record; each task takes about an hour. A task whose .tbl exists is skipped, and
# every task is skipped when a later step's input is supplied (PHMMER_HITS, MAPS, SAMAP_JOINT or SAMAP_FINAL; see
# config.sh). Each task first checks the proteomes against the checksums config.sh pins.
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

# config.sh sets $OUT, the proteomes and their checksums, $PHMMER, exists, pinned, verify_sha256 and skip_if_unused.
source "$HERE/config.sh"

N_CHUNKS=32
HOMOLOGY=$OUT/0_homology

# Sets DIRECTION, QUERY (the proteome split into chunks), TARGET, I (the chunk), CHUNK, TBL and PHMMER_LOG for
# array task $1.
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

# The array task ID, which sbatch sets, is the chunk this task searches (1-64, as above).
TASK=${SLURM_ARRAY_TASK_ID:?submit with sbatch, or set SLURM_ARRAY_TASK_ID to a task, 1-64}
task_paths "$TASK"
if [ -e "$TBL" ]; then
    echo "$TBL exists; skipped"
    exit 0
fi
verify_sha256 "$HUMAN_PROTEOME" "$HUMAN_PROTEOME_SHA256"
verify_sha256 "$MOUSE_PROTEOME" "$MOUSE_PROTEOME_SHA256"
mkdir -p "$(dirname "$CHUNK")" "$(dirname "$TBL")"

# The chunk: awk numbers the proteins at their ">" lines and keeps every line of proteins I, I+32, ...
awk -v n="$N_CHUNKS" -v i="$I" '/^>/ { k++ } k > 0 && (k - 1) % n == i - 1' "$QUERY" > "$CHUNK.partial"
mv "$CHUNK.partial" "$CHUNK"

# The search of the run of record (E-value threshold 1e-5, one thread per core). --noali leaves the alignments out
# of the summary output (about 1 GB per chunk in the run of record); the hit table is the same either way.
"$PHMMER" --cpu "${SLURM_CPUS_PER_TASK:-4}" -E 1e-5 --noali --tblout "$TBL.partial" "$CHUNK" "$TARGET" > "$PHMMER_LOG"
mv "$TBL.partial" "$TBL"
echo "$TBL: $(grep -cv '^#' "$TBL") hits for $(grep -c '^>' "$CHUNK") query proteins"
