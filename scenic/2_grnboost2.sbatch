#!/bin/bash
# Step 2 of 6: one GRNBoost2 run per array task, task n being row n of runs.tsv (1-100: mouse, 101-200: human; seeds
# 1000-1099 in each species).
#   The run's sketch: scripts/scenic_make_loom.py draws 15,000 nuclei from step 1's pooled h5ad by geometric sketching
#     on obsm['X_scVI'], seeded with the run's seed, and writes
#     $OUT/scenic_seeds/<species>/looms_sketch15k/<species>_dsc_sketch15k_seed<seed>.loom (and .manifest.json).
#   GRNBoost2: scripts/grn_only.py runs pySCENIC's arboreto_with_multiprocessing.py on that loom with the run's seed
#     and writes $OUT/scenic_seeds/<species>/adjacencies_sketch15k/adjacencies_seed<seed>.tsv (TF, target, importance).
#   Logs: $OUT/logs/2_sketch_<species>_seed<seed>.log, 2_grnboost2_<species>_seed<seed>.log,
#         logs/2_grnboost2_<job ID>_<task>.out
#   Settings: every setting of each Python call, library defaults included, in $OUT/logs/<log name>.settings.json
#     (scripts/run_with_settings.py; see README.md)
#     For GRNBoost2 this includes arboreto's options, the call made for each target gene and every parameter of
#     the gradient-boosting regressor.
# GRNBoost2 took 4 to 7.5 hours per run on 16 cores in September 2026. A run whose adjacencies exist is skipped, and a
# sketch that exists is reused; both are written under a temporary name and renamed when complete.
#
# Usage:
#   bash 2_grnboost2.sbatch --check                    per run: whether its sketch and adjacencies exist; run nothing
#   sbatch 2_grnboost2.sbatch                          all 200 runs
#   sbatch --array=1,101 2_grnboost2.sbatch            seed 1000 of each species
#   SLURM_ARRAY_TASK_ID=101 bash 2_grnboost2.sbatch    human seed 1000, in the current shell (an interactive session)

#SBATCH --job-name=2_grnboost2
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=16  # GRNBoost2 workers (--n-cpu below)
#SBATCH --mem=64G
#SBATCH --time=12:00:00
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

# config.sh sets $OUT, $PY_SCENIC, $ARBORETO, $SKETCH_N, $TAG, N_RUNS, species_paths <species>, and field <row>
# <column>, which prints one value from runs.tsv.
source "$HERE/config.sh"

# Sets SP, SEED, SP_DIR, TF_LIST, SKETCH (the run's loom) and ADJ (its adjacencies) for row $1 of runs.tsv.
run_paths() {
    species_paths "$(field "$1" species)"
    SEED=$(field "$1" seed)
    SKETCH=$SP_DIR/looms_$TAG/${SP}_dsc_${TAG}_seed$SEED.loom
    ADJ=$SP_DIR/adjacencies_$TAG/adjacencies_seed$SEED.tsv
}

if [ "${1:-}" = "--check" ]; then
    echo "pooled h5ad: mouse $(exists "$OUT/scenic_seeds/mouse/mouse_dsc_pooled.h5ad"), human $(exists "$OUT/scenic_seeds/human/human_dsc_pooled.h5ad")"
    n_done=0
    for row in $(seq 1 "$N_RUNS"); do
        run_paths "$row"
        if [ -s "$ADJ" ]; then
            n_done=$((n_done + 1))
        fi
        echo "$row $SP seed $SEED: sketch $(exists "$SKETCH"), adjacencies $(exists "$ADJ")"
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
if [ -s "$ADJ" ]; then
    echo "$SP seed $SEED: $ADJ exists; nothing to do"
    exit 0
fi
mkdir -p "$(dirname "$SKETCH")" "$OUT/logs"
echo "$SP seed $SEED: started $(date) on $(hostname)"

# The run's sketch, as rule scenic_sketch of the run of record made it: raw counts from the pooled h5ad's
# layers['raw'], geometric sketch on obsm['X_scVI'] seeded with the run's seed, the pooled h5ad left as it is.
if [ ! -e "${SKETCH%.loom}.manifest.json" ]; then
    "$PY_SCENIC" "$HERE/scripts/run_with_settings.py" \
        --settings "$OUT/logs/2_sketch_${SP}_seed$SEED.settings.json" \
        --trace geosketch.sketch:gs_gap \
        -- "$HERE/scripts/scenic_make_loom.py" \
        --h5ad "$SP_DIR/${SP}_dsc_pooled.h5ad" \
        --species "$SP" \
        --tf-list "$TF_LIST" \
        --out-dir "$SP_DIR" \
        --skip-h5ad \
        --raw-layer raw \
        --sketch-n "$SKETCH_N" \
        --sketch-seed "$SEED" \
        --celltype-col cell_subtype \
        --out-loom "$SKETCH" \
        --out-manifest "${SKETCH%.loom}.manifest.json" \
        > "$OUT/logs/2_sketch_${SP}_seed$SEED.log" 2>&1
fi

# GRNBoost2 with the run's seed (rule scenic_grn). The arboreto script is the one the run of record found next to the
# environment's Python; it is named here so that no other copy on the PATH is picked up.
"$PY_SCENIC" "$HERE/scripts/run_with_settings.py" \
    --settings "$OUT/logs/2_grnboost2_${SP}_seed$SEED.settings.json" \
    -- "$HERE/scripts/grn_only.py" \
    --loom "$SKETCH" \
    --tf-list "$TF_LIST" \
    --out "$ADJ" \
    --seed "$SEED" \
    --n-cpu 16 \
    --arboreto-script "$ARBORETO" \
    > "$OUT/logs/2_grnboost2_${SP}_seed$SEED.log" 2>&1

# The log's last line (links written) also goes to the SLURM log.
tail -n 1 "$OUT/logs/2_grnboost2_${SP}_seed$SEED.log"
echo "$SP seed $SEED: finished $(date)"
