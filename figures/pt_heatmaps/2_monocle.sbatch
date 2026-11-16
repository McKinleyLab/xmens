#!/bin/bash
# Step 2 of 3: Monocle 2 itself, from the raw counts to the z matrix of every heatmap, per species.
#   scripts/monocle_z.R builds a CellDataSet from step 1's counts, runs estimateSizeFactors() and
#     estimateDispersions(), smooths every gene drawn in the heatmaps with genSmoothCurves() (the function
#     plot_pseudotime_heatmap() uses) and turns the curves into z with plot_pseudotime_heatmap()'s own lines. Mouse
#     also gives Fig. 4H's z matrix (induced nuclei) and plot_pseudotime_heatmap()'s own drawing of that panel.
#   Reads:  $OUT/inputs/ (step 1); monocle from $MONOCLE_LIB (step 0)
#   Writes: $OUT/z/
#   Logs:   $OUT/logs/2_monocle_<species>.log and logs/2_monocle_<job ID>.out
# A species whose z_<species>.tsv exists is skipped.
#
# Usage:
#   sbatch 2_monocle.sbatch                   both species, one after the other
#   SPECIES=human sbatch 2_monocle.sbatch     one species (submit one per species to run them at the same time)

#SBATCH --job-name=2_monocle
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=180G  # estimateDispersions() densifies genes x nuclei: peaks 99 GB mouse, 127 GB human
#SBATCH --time=02:00:00  # 7-10 min mouse, 9-11 min human
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

# config.sh sets $OUT, $R_MODULE, $MONOCLE_LIB and the cache, library and temporary-file locations.
source "$HERE/config.sh"

# The species still to do; a species whose z_<species>.tsv exists is skipped.
todo=""
for species in ${SPECIES:-mouse human}; do
    if [ -e "$OUT/z/z_$species.tsv" ]; then
        echo "$species: $OUT/z/z_$species.tsv exists; skipped"
    else
        todo="$todo $species"
    fi
done
if [ -z "$todo" ]; then
    exit 0
fi

# R 4.4.1 from the cluster's modules (Lmod does not run under set -u); monocle's library is searched first.
set +u
source /etc/profile.d/lmod.sh
module load "$R_MODULE"
set -u
export R_LIBS=$MONOCLE_LIB

mkdir -p "$OUT/z" "$OUT/logs"
for species in $todo; do
    Rscript "$HERE/scripts/monocle_z.R" "$species" "$OUT/inputs" "$OUT/z" \
        2>&1 | tee "$OUT/logs/2_monocle_$species.log"
done
