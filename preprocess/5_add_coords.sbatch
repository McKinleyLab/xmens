#!/bin/bash
# Step 5 of 5: add Andy's Slide-tags positions to the QC'd objects, next to Nick's.
#   Reads:  $OUT/qc_passed_{mouse,human}.h5ad (step 4), Nick's CSVs ($COORDS_DIR) and Andy's positions ($ANDY_DIR)
#   Writes: $OUT/with_coords/qc_passed_{mouse,human}.h5ad, copies with obs x_um_dbscan_nick / y_um_dbscan_nick
#           (Nick's, as in obsm['spatial']) and x_um_dbscan_andy / y_um_dbscan_andy (Andy's, also
#           obsm['spatial_andy']), and $OUT/with_coords/coordinate_summary.tsv
#   Log:    logs/5_add_coords_<job ID>.out
# Stops before copying anything if with_coords/ exists or if the objects' positions differ from Nick's CSVs.
#
# Usage:
#   sbatch 5_add_coords.sbatch

#SBATCH --job-name=5_add_coords
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out  # job name and job ID; relative to the folder sbatch was run from

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# The pipeline folder holds config.sh and samples.tsv. Run with bash, this script sits in it; under sbatch, SLURM
# runs a copy from its spool folder, and the pipeline folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the preprocess folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, $PY, $COORDS_DIR and $ANDY_DIR.
source "$HERE/config.sh"

"$PY" "$HERE/scripts/add_coords.py" "$OUT" "$OUT/with_coords" \
    --nick-dir "$COORDS_DIR" \
    --andy-dir "$ANDY_DIR"
