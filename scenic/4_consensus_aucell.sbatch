#!/bin/bash
# Step 4 of 6: the 100 runs -> consensus regulons and their activity in every nucleus. One array task per species
# (1 mouse, 2 human).
#   Consensus: scripts/build_consensus_aucell.py reads the regulons of the 100 runs from $REGULONS_FROM (config.sh; the
#     run folder unless set). It keeps a TF that formed a regulon in at least 80 of the 100 runs and has at least 10
#     targets, and a target whose link to the TF recurs in at least <floor> runs, weighted by the fraction of the TF's
#     runs that contain it. It runs twice, each into its own folder of $OUT/scenic_seeds/<species>/:
#       consensus_sketch15k_weighted_soft_occ6/  floor 6: the regulons of record (193 mouse, 141 human)
#       consensus_sketch15k_weighted_soft/       floor 2: Fig. S28 (step 6) reads its links
#     Each holds consensus_regulons.pkl, recurrence_report.tsv, recurrence_summary.txt, recurrence_hist.png, and
#     auc_mtx.csv and auc_by_timepoint.tsv: the script's own AUCell over all genes of the pooled h5ad, kept as a
#     comparator. The activity of record is the next part.
#   AUCell of record: scripts/scenic_regulons_to_yaml.py writes the floor-6 regulons as YAML (pyscenic aucell cannot read
#     the pickle) and checks that they reload unchanged; pySCENIC's command-line AUCell then scores them on step 1's
#     pooled loom (genes detected in at least 1% of nuclei), ranking seed 0, AUC threshold 0.05, into
#     consensus_sketch15k_weighted_soft_occ6_pyscenic_cli/:
#       auc_mtx_weighted.csv    --weights (targets weighted by recurrence): the activity of record
#       auc_mtx_unweighted.csv  without --weights, for comparison
#   Logs: $OUT/logs/4_consensus_<species>_floor6.log, _floor2.log, 4_aucell_<species>.log,
#         logs/4_consensus_aucell_<job ID>_<task>.out
#   Settings: every setting of each Python call, library defaults included, in $OUT/logs/<log name>.settings.json
#     (scripts/run_with_settings.py; see README.md)
# Each folder is written as <folder>.partial.<job ID> and renamed when complete; a folder that exists is skipped. A
# failed task leaves its .partial folder behind, to be looked at and removed by hand.
#
# Usage:
#   bash 4_consensus_aucell.sbatch --check                  per species: inputs and outputs found; run nothing
#   sbatch 4_consensus_aucell.sbatch                        both species
#   sbatch --array=2 4_consensus_aucell.sbatch              human only
#   SLURM_ARRAY_TASK_ID=2 bash 4_consensus_aucell.sbatch    human, in the current shell (an interactive session)

#SBATCH --job-name=4_consensus_aucell
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=8  # AUCell workers (--n-cpu and --num_workers below)
#SBATCH --mem=128G  # the consensus script's AUCell holds the dense nuclei x genes matrix (43,197 x 33,698 in mouse)
#SBATCH --time=08:00:00
#SBATCH --array=1-2  # 1 mouse, 2 human
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

# config.sh sets $OUT, $REGULONS_FROM, $PY_SCENIC, $PYSCENIC_CLI, $SKETCH, $TAG, $MIN_OCCURRENCE, the folder names
# $CONS6, $CONS2 and $AUCELL, and species_paths <species>.
source "$HERE/config.sh"

if [ "${1:-}" = "--check" ]; then
    echo "regulons of the 100 runs from: $REGULONS_FROM"
    for task in 1 2; do
        species_paths "$(task_species "$task")"
        n_pkl=$(ls "$REGULONS_FROM/scenic_seeds/$SP/regulons_$TAG"/regulons_seed*.pkl 2>/dev/null | wc -l | tr -d ' ' || true)
        echo "$task $SP: $n_pkl run regulon files, pooled h5ad $(exists "$SP_DIR/${SP}_dsc_pooled.h5ad")," \
            "pooled loom $(exists "$SP_DIR/${SP}_dsc_pooled.loom")"
        echo "    floor 6 $(exists "$SP_DIR/$CONS6"), floor 2 $(exists "$SP_DIR/$CONS2"), AUCell $(exists "$SP_DIR/$AUCELL")"
    done
    exit 0
fi

TASK=${SLURM_ARRAY_TASK_ID:?submit with sbatch, or set SLURM_ARRAY_TASK_ID to 1 (mouse) or 2 (human)}
SPECIES=$(task_species "$TASK")  # on its own line, so that a task other than 1 or 2 stops the script here
species_paths "$SPECIES"
JOB=${SLURM_JOB_ID:-$$}
mkdir -p "$OUT/logs"
echo "$SP: started $(date) on $(hostname); regulons of the runs from $REGULONS_FROM"

# consensus <floor> <folder>: the command of the run of record (rule scenic_consensus; the floor-6 folder was made by
# the same command on 14 Oct 2026), written into <folder>.partial.<job ID> and renamed when complete.
consensus() {
    local floor=$1 dest=$SP_DIR/$2
    local tmp=$dest.partial.$JOB
    if [ -e "$dest" ]; then
        echo "$SP: $dest exists; skipped"
        return 0
    fi
    "$PY_SCENIC" "$HERE/scripts/run_with_settings.py" \
        --settings "$OUT/logs/4_consensus_${SP}_floor$floor.settings.json" \
        --trace pyscenic.aucell:aucell \
        -- "$HERE/scripts/build_consensus_aucell.py" \
        --species "$SP" \
        --sketch "$SKETCH" \
        --mode weighted_soft \
        --rev "$REGULONS_FROM" \
        --h5ad "$SP_DIR/${SP}_dsc_pooled.h5ad" \
        --timepoint-col timepoint \
        --min-occurrence "$floor" \
        --out-dir "$tmp" \
        --n-cpu 8 \
        > "$OUT/logs/4_consensus_${SP}_floor$floor.log" 2>&1
    mv "$tmp" "$dest"
    grep -E "^kept regulons|^consensus regulons" "$OUT/logs/4_consensus_${SP}_floor$floor.log"
}
consensus "$MIN_OCCURRENCE" "$CONS6"
consensus 2 "$CONS2"

# The AUCell of record (the 26 Oct 2026 run, cliauc_run.sh; rule scenic_aucell_cli): YAML export with a reload check,
# then pyscenic aucell on the pooled loom without and with --weights, ranking seed 0.
DEST=$SP_DIR/$AUCELL
if [ -e "$DEST" ]; then
    echo "$SP: $DEST exists; skipped"
else
    TMP=$DEST.partial.$JOB
    LOG=$OUT/logs/4_aucell_$SP.log
    mkdir -p "$TMP"
    "$PY_SCENIC" "$HERE/scripts/run_with_settings.py" \
        --settings "$OUT/logs/4_aucell_${SP}_yaml.settings.json" \
        -- "$HERE/scripts/scenic_regulons_to_yaml.py" \
        --pkl "$SP_DIR/$CONS6/consensus_regulons.pkl" \
        --yaml "$TMP/consensus_regulons.yaml" \
        > "$LOG" 2>&1
    "$PY_SCENIC" "$HERE/scripts/run_with_settings.py" \
        --settings "$OUT/logs/4_aucell_${SP}_unweighted.settings.json" \
        --trace pyscenic.cli.utils:load_exp_matrix \
        --trace pyscenic.aucell:aucell \
        -- "$PYSCENIC_CLI" aucell "$SP_DIR/${SP}_dsc_pooled.loom" "$TMP/consensus_regulons.yaml" \
        --output "$TMP/auc_mtx_unweighted.csv" \
        --num_workers 8 \
        --seed 0 \
        >> "$LOG" 2>&1
    "$PY_SCENIC" "$HERE/scripts/run_with_settings.py" \
        --settings "$OUT/logs/4_aucell_${SP}_weighted.settings.json" \
        --trace pyscenic.cli.utils:load_exp_matrix \
        --trace pyscenic.aucell:aucell \
        -- "$PYSCENIC_CLI" aucell "$SP_DIR/${SP}_dsc_pooled.loom" "$TMP/consensus_regulons.yaml" \
        --output "$TMP/auc_mtx_weighted.csv" \
        --num_workers 8 \
        --seed 0 \
        --weights \
        >> "$LOG" 2>&1
    mv "$TMP" "$DEST"
    head -n 1 "$LOG"
fi
echo "$SP: finished $(date)"
