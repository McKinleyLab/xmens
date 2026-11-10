#!/bin/bash
# Step 3 of 6: cisTarget pruning, one run per array task, task n being row n of runs.tsv (as in step 2).
#   scripts/ctx_only.py reads the run's adjacencies and sketch loom from $GRN_FROM (config.sh; the run folder unless
#   set), builds pySCENIC's co-expression modules (positive correlations only, dropouts not masked), keeps the modules
#   enriched for a motif annotated to their TF (prune2df with the Python API's defaults: rank threshold 1,500, AUC
#   threshold 0.05, NES 3; mm10 or hg38 rankings, 500 bp and 10 kb, in $DB_DIR) and writes the run's regulons to
#   $OUT/scenic_seeds/<species>/regulons_sketch15k/regulons_seed<seed>.pkl
#   Logs: $OUT/logs/3_cistarget_<species>_seed<seed>.log, logs/3_cistarget_<job ID>_<task>.out
#   Settings: every setting of each Python call, library defaults included, in $OUT/logs/<log name>.settings.json
#     (scripts/run_with_settings.py; see README.md)
# About 5 minutes per run in September 2026. A run whose regulons exist is skipped; the pickle is written under a
# temporary name and renamed when complete.
#
# Usage:
#   bash 3_cistarget.sbatch --check                    per run: whether its inputs and regulons exist; run nothing
#   sbatch 3_cistarget.sbatch                          all 200 runs
#   sbatch --array=1,101 3_cistarget.sbatch            seed 1000 of each species
#   SLURM_ARRAY_TASK_ID=101 bash 3_cistarget.sbatch    human seed 1000, in the current shell (an interactive session)

#SBATCH --job-name=3_cistarget
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=8  # pruning workers (--n-cpu below)
#SBATCH --mem=96G  # the gene rankings are loaded whole
#SBATCH --time=02:00:00
#SBATCH --array=1-200  # one task per row of runs.tsv
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506  # GPU nodes
#SBATCH --output=logs/%x_%A_%a.out  # job name, job ID, task; relative to the folder sbatch was run from

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# The pipeline folder holds config.sh and runs.tsv. Run with bash, this script sits in it; under sbatch, SLURM runs a
# copy from its spool folder, and the pipeline folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the scenic folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, $GRN_FROM, $DB_DIR, $PY_SCENIC, $TAG, N_RUNS, species_paths <species>, and field <row> <column>.
source "$HERE/config.sh"

# Sets SP, SEED, SP_DIR, ADJ and SKETCH (this run's inputs, under $GRN_FROM) and PKL (its output) for row $1 of runs.tsv.
run_paths() {
    species_paths "$(field "$1" species)"
    SEED=$(field "$1" seed)
    ADJ=$GRN_FROM/scenic_seeds/$SP/adjacencies_$TAG/adjacencies_seed$SEED.tsv
    SKETCH=$GRN_FROM/scenic_seeds/$SP/looms_$TAG/${SP}_dsc_${TAG}_seed$SEED.loom
    PKL=$SP_DIR/regulons_$TAG/regulons_seed$SEED.pkl
}

if [ "${1:-}" = "--check" ]; then
    echo "inputs from: $GRN_FROM"
    n_done=0
    for row in $(seq 1 "$N_RUNS"); do
        run_paths "$row"
        if [ -s "$PKL" ]; then
            n_done=$((n_done + 1))
        fi
        echo "$row $SP seed $SEED: adjacencies $(exists "$ADJ"), sketch $(exists "$SKETCH"), regulons $(exists "$PKL")"
    done
    echo "$n_done of $N_RUNS runs done"
    exit 0
fi

# The array task ID, which sbatch sets, is the row of runs.tsv this task runs.
ROW=${SLURM_ARRAY_TASK_ID:?submit with sbatch, or set SLURM_ARRAY_TASK_ID to a row of runs.tsv}
if [ "$ROW" -lt 1 ] || [ "$ROW" -gt "$N_RUNS" ]; then
    echo "runs.tsv has rows 1-$N_RUNS, not $ROW" >&2
    exit 2
fi
run_paths "$ROW"
if [ -s "$PKL" ]; then
    echo "$SP seed $SEED: $PKL exists; nothing to do"
    exit 0
fi
mkdir -p "$OUT/logs"
echo "$SP seed $SEED: started $(date) on $(hostname); inputs from $GRN_FROM"

# The command of the run of record (rule scenic_ctx).
"$PY_SCENIC" "$HERE/scripts/run_with_settings.py" \
    --settings "$OUT/logs/3_cistarget_${SP}_seed$SEED.settings.json" \
    --trace pyscenic.utils:modules_from_adjacencies \
    --trace pyscenic.prune:prune2df \
    --trace pyscenic.prune:df2regulons \
    -- "$HERE/scripts/ctx_only.py" \
    --adjacencies "$ADJ" \
    --loom "$SKETCH" \
    --db-dir "$DB_DIR" \
    --species "$SP" \
    --out "$PKL" \
    --n-cpu 8 \
    > "$OUT/logs/3_cistarget_${SP}_seed$SEED.log" 2>&1

# The number of regulons written also goes to the SLURM log.
grep "Wrote" "$OUT/logs/3_cistarget_${SP}_seed$SEED.log"
echo "$SP seed $SEED: finished $(date)"
