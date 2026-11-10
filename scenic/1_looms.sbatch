#!/bin/bash
# Step 1 of 6: fibroblast/decidual object -> the pooled objects every later step reads. One array task per species
# (1 mouse, 2 human).
#   scripts/scenic_make_loom.py reads $SCENIC_INPUT_MOUSE or $SCENIC_INPUT_HUMAN (config.sh), takes the raw counts
#   (layers['counts']) and writes, in $OUT/scenic_seeds/<species>/:
#     <species>_dsc_pooled.h5ad           every nucleus and gene, raw counts (X and layers['raw']); step 2 sketches it
#                                         and step 4's consensus script scores it
#     <species>_dsc_pooled.loom           genes detected in at least 1% of nuclei (14,245 mouse, 15,011 human); step 4's
#                                         AUCell of record scores it
#     <species>_dsc_pooled.manifest.json  numbers of nuclei, genes and TFs, and the input's sha256
#   Logs: $OUT/logs/1_looms_<species>.log, logs/1_looms_<job ID>_<task>.out
#   Settings: every setting of each Python call, library defaults included, in $OUT/logs/<log name>.settings.json
#     (scripts/run_with_settings.py; see README.md)
# A species whose manifest exists is skipped. The script writes each file under a temporary name, the manifest last.
#
# Usage:
#   bash 1_looms.sbatch --check                  print each species' input and whether it is done; run nothing
#   sbatch 1_looms.sbatch                        both species
#   sbatch --array=2 1_looms.sbatch              human only
#   SLURM_ARRAY_TASK_ID=2 bash 1_looms.sbatch    human, in the current shell (an interactive session)

#SBATCH --job-name=1_looms
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G  # the input object, 1.5 GB on disk, is read whole
#SBATCH --time=02:00:00
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

# config.sh sets $OUT, $PY_SCENIC, the inputs, and species_paths <species>, which sets SP, SP_DIR, TF_LIST and
# SCENIC_INPUT.
source "$HERE/config.sh"

if [ "${1:-}" = "--check" ]; then
    for task in 1 2; do
        species_paths "$(task_species "$task")"
        state="to run"
        if [ -e "$SP_DIR/${SP}_dsc_pooled.manifest.json" ]; then
            state="done"
        fi
        echo "$task $SP: input $(exists "$SCENIC_INPUT"), TF list $(exists "$TF_LIST"), $state"
        echo "    input: $SCENIC_INPUT"
    done
    exit 0
fi

TASK=${SLURM_ARRAY_TASK_ID:?submit with sbatch, or set SLURM_ARRAY_TASK_ID to 1 (mouse) or 2 (human)}
SPECIES=$(task_species "$TASK")  # on its own line, so that a task other than 1 or 2 stops the script here
species_paths "$SPECIES"
MANIFEST=$SP_DIR/${SP}_dsc_pooled.manifest.json
if [ -e "$MANIFEST" ]; then
    echo "$SP: $MANIFEST exists; nothing to do"
    exit 0
fi
mkdir -p "$SP_DIR" "$OUT/logs"

# The command of the run of record (rule scenic_pooled, xmens_final/workflow/rules/24_scenic.smk): raw counts from
# layers['counts'], no sketch. --celltype-col only labels the composition table in the log.
echo "$SP: started $(date) on $(hostname)"
"$PY_SCENIC" "$HERE/scripts/run_with_settings.py" \
    --settings "$OUT/logs/1_looms_$SP.settings.json" \
    -- "$HERE/scripts/scenic_make_loom.py" \
    --h5ad "$SCENIC_INPUT" \
    --species "$SP" \
    --tf-list "$TF_LIST" \
    --out-dir "$SP_DIR" \
    --raw-layer "$RAW_LAYER" \
    --sketch-n 0 \
    --celltype-col cell_subtype \
    --out-loom "$SP_DIR/${SP}_dsc_pooled.loom" \
    --out-h5ad "$SP_DIR/${SP}_dsc_pooled.h5ad" \
    --out-manifest "$MANIFEST" \
    > "$OUT/logs/1_looms_$SP.log" 2>&1

# The gene filter and the files written also go to the SLURM log.
grep -E "^genes:|^wrote" "$OUT/logs/1_looms_$SP.log"
echo "$SP: finished $(date)"
