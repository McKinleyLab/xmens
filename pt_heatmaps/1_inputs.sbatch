#!/bin/bash
# Step 1 of 3: final objects, MI tables and gene lists -> the inputs of Monocle (step 2) and of the figures (step 3).
#   scripts/heatmap_inputs.py writes, per species, the raw counts of every gene for the nuclei with a pseudotime
#     (42,290 mouse, 47,533 human) and the table of those nuclei, plus the genes of each heatmap: the Yang et al. (2023)
#     panel (Fig. 4H), the 28 genes of Fig. 5F, and the top MI genes of each species (big MI heatmaps).
#   Reads:  $H5AD_MOUSE, $H5AD_HUMAN, $MI_MOUSE, $MI_HUMAN, $PAIRS, $FIG5F_GENES (config.sh) and
#           resources/yang2023_panel.tsv
#   Writes: $OUT/inputs/
#   Logs:   $OUT/logs/1_inputs.log and logs/1_inputs_<job ID>.out
# Skipped if $OUT/inputs exists.
#
# Usage:
#   bash 1_inputs.sbatch --check    print the inputs and whether they exist; run nothing
#   sbatch 1_inputs.sbatch

#SBATCH --job-name=1_inputs
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
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
    echo "submit this from the xmens_pt_heatmaps folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, $PY, the input paths and exists <path>.
source "$HERE/config.sh"

if [ "${1:-}" = "--check" ]; then
    for path in "$H5AD_MOUSE" "$H5AD_HUMAN" "$MI_MOUSE" "$MI_HUMAN" "$PAIRS" "$FIG5F_GENES" \
            "$HERE/resources/yang2023_panel.tsv"; do
        echo "$(exists "$path")  $path"
    done
    echo "outputs $(exists "$OUT/inputs"): $OUT/inputs"
    exit 0
fi

if [ -e "$OUT/inputs" ]; then
    echo "$OUT/inputs exists; skipped"
    exit 0
fi
mkdir -p "$OUT/logs"
"$PY" "$HERE/scripts/heatmap_inputs.py" \
    --h5ad-mouse "$H5AD_MOUSE" \
    --h5ad-human "$H5AD_HUMAN" \
    --mi-mouse "$MI_MOUSE" \
    --mi-human "$MI_HUMAN" \
    --pairs "$PAIRS" \
    --fig5f-genes "$FIG5F_GENES" \
    --yang "$HERE/resources/yang2023_panel.tsv" \
    --out-dir "$OUT/inputs" \
    2>&1 | tee "$OUT/logs/1_inputs.log"
