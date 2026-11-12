#!/bin/bash
# Step 6 of 6: the 28 genes of the cross-species pseudotime heatmap (Fig. 5F): from the 273 one-to-one orthologs in the
# top 10% by MI in both species, the 14 early and 14 late genes whose human and mouse pseudotime profiles agree best.
#   Reads:  $OUT/mi_mouse.csv and $OUT/mi_human.csv (step 3), $OUT/overlap/ortholog_pairs_one2one.csv (step 4),
#           $H5AD_MOUSE and $H5AD_HUMAN (config.sh)
#   Writes: $OUT/fig5f/conserved_profiles_selected_genes.csv (the 28 genes, in drawing order),
#           conserved_profiles_pool_concordance.csv, selection_summary.json and percell_<species>_pool.tsv.gz
#   Logs:   $OUT/logs/6_fig5f_genes.log and logs/6_fig5f_genes_<job ID>.out
# Skipped if $OUT/fig5f exists. The heatmap itself is drawn from these files by the Fig. 5F figure scripts.
#
# Usage:
#   sbatch 6_fig5f_genes.sbatch

#SBATCH --job-name=6_fig5f_genes
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G  # one object's counts layer at a time; the test run took 3 min and 0.5 GB
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

# config.sh sets $OUT, $PY, $H5AD_MOUSE and $H5AD_HUMAN.
source "$HERE/config.sh"

if [ -e "$OUT/fig5f" ]; then
    echo "$OUT/fig5f exists; skipped"
    exit 0
fi
mkdir -p "$OUT/logs"
"$PY" "$HERE/scripts/fig5f_genes.py" \
    --mi-mouse "$OUT/mi_mouse.csv" \
    --mi-human "$OUT/mi_human.csv" \
    --pairs "$OUT/overlap/ortholog_pairs_one2one.csv" \
    --h5ad-mouse "$H5AD_MOUSE" \
    --h5ad-human "$H5AD_HUMAN" \
    --counts-layer "$COUNTS_LAYER" \
    --out-dir "$OUT/fig5f" \
    2>&1 | tee "$OUT/logs/6_fig5f_genes.log"
