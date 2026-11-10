#!/bin/bash
# Step 3 of 3: the figures, by running the notebook xmens_pt_heatmaps.ipynb (its first cell finds steps 1-2 done; then
# one cell per figure). Colours and layout only: Monocle's z (step 2), clipped at +/-3, through RdBu_r.
#   Fig. 4H   the Yang et al. (2023) panel, rows in published order: the 30 x 59.887 mm box, the 5 pt colorbar and
#             the labelled panel
#   Fig. 5F   the 28 conserved genes, human (left) and mouse (right)
#   MI        every top MI gene of each species, rows in order of pseudotime centroid, with the gene-order tables
#   Reads:  $OUT/inputs/ (step 1) and $OUT/z/ (step 2)
#   Writes: $OUT/figures/fig4h/, fig5f/ and mi/, and $OUT/xmens_pt_heatmaps.executed.ipynb (the notebook as run)
#   Log:    logs/3_figures_<job ID>.out
# The figures are drawn again each time this runs.
#
# Usage:
#   sbatch 3_figures.sbatch

#SBATCH --job-name=3_figures
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
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
    echo "submit this from the xmens_pt_heatmaps folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT and $PY, and keeps the notebook kernel's files out of the home folder.
source "$HERE/config.sh"

# The notebook would start Monocle itself if step 2 were not done, in this job's 8 GB.
for species in mouse human; do
    if [ ! -e "$OUT/z/z_$species.tsv" ]; then
        echo "no $OUT/z/z_$species.tsv: run step 2 first" >&2
        exit 1
    fi
done

# The rsc_uv environment's own kernel (python3, the Python running Jupyter here), so no kernel needs installing.
cd "$HERE"
"$PY" -m jupyter nbconvert --to notebook --execute xmens_pt_heatmaps.ipynb \
    --ExecutePreprocessor.kernel_name=python3 --ExecutePreprocessor.timeout=1800 \
    --output "$OUT/xmens_pt_heatmaps.executed.ipynb"
ls -l "$OUT/figures"/*
