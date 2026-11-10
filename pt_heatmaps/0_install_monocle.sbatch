#!/bin/bash
# Step 0 (once): install Monocle 2 for R 4.4.1.
#   scripts/install_monocle.R installs monocle 2.34.0 (Bioconductor 3.20) and the dependencies the lab's R library lacks
#     into $MONOCLE_LIB (config.sh), plus dplyr 1.1.4: Monocle's estimateDispersions() calls functions that dplyr 1.2
#     (the lab's library has 1.2.1) made defunct. Packages already in the lab's library are used as they are; nothing
#     there is installed, updated or removed. Ends by running Monocle's dispersion estimate on simulated data.
#   Writes: $MONOCLE_LIB/ and $MONOCLE_LIB/monocle_install_sessionInfo.txt
#   Log:    logs/0_install_monocle_<job ID>.out
# Downloads from Bioconductor and CRAN, and compiles a few packages. Installs only what is missing.
#
# Usage:
#   sbatch 0_install_monocle.sbatch

#SBATCH --job-name=0_install_monocle
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
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

# config.sh sets $MONOCLE_LIB, $R_MODULE and the cache, library and temporary-file locations. It needs $OUT, which
# this step does not use.
OUT=${OUT:-unused}
source "$HERE/config.sh"

# R 4.4.1 from the cluster's modules (Lmod does not run under set -u).
set +u
source /etc/profile.d/lmod.sh
module load "$R_MODULE"
set -u

Rscript "$HERE/scripts/install_monocle.R" "$MONOCLE_LIB"
# The check loads monocle as step 2 does: its library first on the search path.
R_LIBS=$MONOCLE_LIB Rscript "$HERE/scripts/check_monocle.R" "$MONOCLE_LIB"
