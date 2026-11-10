#!/bin/bash
# Step 2 of 5: the SAMap fit.
#   notebooks/2_samap.ipynb (kernel samap) builds each species' SAM with SAMap's defaults and runs SAMap (sc-samap
#   3.0.1) for three iterations with the hmmer homology maps: arm vanilla_raw of the run of record.
#   Reads:  $SAMAP_INPUTS/{hs,mm}_fib_dec_all_raw.h5ad (step 1), $MAPS/hsmm/{hs_to_mm,mm_to_hs}.txt (step 0b)
#   Writes: $OUT/2_samap/ samap_joint.h5ad (graph, UMAP, SAM-processed expression), samap.pkl (the SAMAP object,
#           about 11 GB), eval.json, mapping_{celltype,subtype,timepoint}.csv, figures/, meta.json
#   Logs:   $OUT/notebooks/2_samap.ipynb and logs/2_samap_<job ID>.out
# Resources and run time (about 50 min) are those of the run of record. Skipped when $SAMAP_JOINT exists (by default
# $OUT/2_samap/samap_joint.h5ad, which appears when the notebook renames 2_samap.partial at the end; then it must have
# been made from $SAMAP_INPUTS and $MAPS), or when $SAMAP_FINAL, the next step's output, is supplied (see config.sh).
#
# Usage:
#   bash 2_samap.sbatch --check
#   sbatch 2_samap.sbatch

#SBATCH --job-name=2_samap
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=8  # SAMap's ncpus
#SBATCH --mem=128G
#SBATCH --time=06:00:00
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
    echo "submit this from the samap folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, $SAMAP_INPUTS, $SAMAP_JOINT, $MAPS, run_notebook and exists.
source "$HERE/config.sh"

if [ "${1:-}" = "--check" ]; then
    for f in "$SAMAP_INPUTS/hs_fib_dec_all_raw.h5ad" "$SAMAP_INPUTS/mm_fib_dec_all_raw.h5ad" \
             "$MAPS/hsmm/hs_to_mm.txt" "$MAPS/hsmm/mm_to_hs.txt"; do
        echo "reads:  $f ($(exists "$f"))"
    done
    echo "writes: $SAMAP_JOINT ($(exists "$SAMAP_JOINT")) and the rest of $OUT/2_samap/"
    exit 0
fi

skip_if_unused 2 SAMAP_FINAL
if [ -e "$SAMAP_JOINT" ]; then
    if [ "$SAMAP_JOINT" = "$DEFAULT_SAMAP_JOINT" ]; then
        made_from "$OUT/2_samap/meta.json" inputs/hs/path "$SAMAP_INPUTS/hs_fib_dec_all_raw.h5ad"
        made_from "$OUT/2_samap/meta.json" maps/hs_to_mm.txt/path "$MAPS/hsmm/hs_to_mm.txt"
    fi
    echo "$SAMAP_JOINT exists; skipped"
    exit 0
fi
if [ "$SAMAP_JOINT" != "$OUT/2_samap/samap_joint.h5ad" ]; then
    echo "SAMAP_JOINT is set to $SAMAP_JOINT, which does not exist; step 2 writes only to $OUT/2_samap" >&2
    exit 1
fi

run_notebook "$HERE/notebooks/2_samap.ipynb" samap "$OUT/notebooks/2_samap.ipynb" \
    -p OUT "$OUT" \
    -p INPUTS "$SAMAP_INPUTS" \
    -p MAPS "$MAPS" \
    -p NCPUS "${SLURM_CPUS_PER_TASK:-8}"
