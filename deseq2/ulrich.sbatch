#!/bin/bash
# The Ulrich et al. 2024 atlas fit: decidual stroma (dS1-3) vs C7-Fibroblast, within donor, all 10 donors with both.
#   ulrich_decidual_vs_C7_fibroblast.ipynb (kernel rsc): one pseudobulk per donor and state from the raw counts of the
#   CELLxGENE object, PyDESeq2 Wald test under ~ donor + state, checked against the table of record.
#   Reads:  $ULRICH_H5AD, $RESOURCES/ulrich_feature_manifest.tsv, and in $RECORD
#           ulrich_dS_vs_C7fib_all_donors_all_genes.csv.gz and ulrich_dS_vs_C7fib_all_donors_metadata.tsv
#   Writes: $OUT/ulrich/ ulrich_dS_vs_C7fib_all_donors_all_genes.csv, ulrich_dS_vs_C7fib_all_donors_metadata.tsv
#           (pseudobulks) and ulrich_dS_vs_C7fib_all_donors_counts.tsv (pseudobulk counts)
#   Logs:   $OUT/notebooks/ulrich_decidual_vs_C7_fibroblast.ipynb and logs/ulrich_<job ID>.out
# Skipped when $OUT/ulrich exists.
#
# Usage:
#   bash ulrich.sbatch --check
#   sbatch ulrich.sbatch

#SBATCH --job-name=ulrich
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

STEP_OUT=$OUT/ulrich
MANIFEST=$RESOURCES/ulrich_feature_manifest.tsv
RECORD_TABLE=$RECORD/ulrich_dS_vs_C7fib_all_donors_all_genes.csv.gz
RECORD_PSEUDOBULKS=$RECORD/ulrich_dS_vs_C7fib_all_donors_metadata.tsv
if [ "${1:-}" = "--check" ]; then
    for f in "$ULRICH_H5AD" "$MANIFEST" "$RECORD_TABLE" "$RECORD_PSEUDOBULKS"; do
        echo "reads:  $f ($(exists "$f"))"
    done
    echo "writes: $STEP_OUT ($(exists "$STEP_OUT"))"
    exit 0
fi

if [ -e "$STEP_OUT" ]; then
    echo "$STEP_OUT exists; skipped"
    exit 0
fi

run_notebook "$HERE/ulrich_decidual_vs_C7_fibroblast.ipynb" rsc "$OUT/notebooks/ulrich_decidual_vs_C7_fibroblast.ipynb" \
    -p OUT_DIR "$STEP_OUT" \
    -p N_CPUS "${SLURM_CPUS_PER_TASK:-8}" \
    -p H5AD "$ULRICH_H5AD" \
    -p FEATURE_MANIFEST "$MANIFEST" \
    -p RECORD "$RECORD_TABLE" \
    -p RECORD_PSEUDOBULKS "$RECORD_PSEUDOBULKS"
