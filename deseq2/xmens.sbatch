#!/bin/bash
# X-Mens decidual vs fibroblast nuclei: the pseudobulk DESeq2 fit of record of 27 October 2026.
#   xmens_decidual_vs_fibroblast.ipynb (kernel rsc): sums the decidual and fibroblast nuclei of final_mouse.h5ad into
#   11 pseudobulks (both cell types at the five induced timepoints, and the uninduced fibroblasts), fits PyDESeq2
#   ~ timepoint + group, and checks the results table and the pseudobulks against the record.
#   Reads:  $XMENS_H5AD, $RESOURCES/xmens_feature_manifest.tsv, $RECORD/xmens_broad_*_with_uninduced.*
#   Writes: $OUT/xmens/ xmens_broad_all_genes_with_uninduced.csv, xmens_broad_carrier_counts_with_uninduced.tsv,
#           xmens_broad_metadata_with_uninduced.tsv, xmens_broad_genes.tsv, fit_diagnostics_with_uninduced.json
#   Logs:   $OUT/notebooks/xmens_decidual_vs_fibroblast.ipynb and logs/xmens_<job ID>.out
# Skipped when $OUT/xmens exists.
#
# Usage:
#   bash xmens.sbatch --check
#   sbatch xmens.sbatch

#SBATCH --job-name=xmens
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out  # job name and job ID; relative to the folder sbatch was run from

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# Run with bash, this script sits in the notebook folder; under sbatch, SLURM runs a copy from its spool folder, and
# the notebook folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the deseq2 notebook folder (no config.sh in $HERE)" >&2
    exit 2
fi
source "$HERE/config.sh"

STEP_OUT=$OUT/xmens
MANIFEST=$RESOURCES/xmens_feature_manifest.tsv
RECORD_TABLE=$RECORD/xmens_broad_all_genes_with_uninduced.csv.gz
RECORD_PSEUDOBULKS=$RECORD/xmens_broad_metadata_with_uninduced.tsv
RECORD_COUNTS=$RECORD/xmens_broad_carrier_counts_with_uninduced.tsv.gz
if [ "${1:-}" = "--check" ]; then
    for f in "$XMENS_H5AD" "$MANIFEST" "$RECORD_TABLE" "$RECORD_PSEUDOBULKS" "$RECORD_COUNTS"; do
        echo "reads:  $f ($(exists "$f"))"
    done
    echo "writes: $STEP_OUT ($(exists "$STEP_OUT"))"
    exit 0
fi

if [ -e "$STEP_OUT" ]; then
    echo "$STEP_OUT exists; skipped"
    exit 0
fi

run_notebook "$HERE/xmens_decidual_vs_fibroblast.ipynb" rsc "$OUT/notebooks/xmens_decidual_vs_fibroblast.ipynb" \
    -p OUT_DIR "$STEP_OUT" \
    -p N_CPUS "${SLURM_CPUS_PER_TASK:-8}" \
    -p H5AD "$XMENS_H5AD" \
    -p FEATURE_MANIFEST "$MANIFEST" \
    -p RECORD "$RECORD_TABLE" \
    -p RECORD_PSEUDOBULKS "$RECORD_PSEUDOBULKS" \
    -p RECORD_COUNTS "$RECORD_COUNTS"
