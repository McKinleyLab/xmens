#!/bin/bash
# Human endometrium (Slide-tags, 2 donors): decidual vs fibroblast nuclei, the pseudobulk DESeq2 fit of record.
#   human_decidual_vs_fibroblast.ipynb (kernel rsc): sums layers/counts into one pseudobulk per donor and arm (4),
#   tests decidual vs fibroblast with PyDESeq2 (~ donor + group) and checks the table against the record.
#   Reads:  $HUMAN_H5AD, $RECORD/human_decidual_vs_fibroblast_all_genes.csv.gz
#   Writes: $OUT/human/ human_decidual_vs_fibroblast_all_genes.csv, pseudobulks.csv, pseudobulk_counts.csv.gz,
#           pseudobulk_detected.csv.gz, check_against_record.csv
#   Logs:   $OUT/notebooks/human_decidual_vs_fibroblast.ipynb and logs/human_<job ID>.out
# Skipped when $OUT/human exists.
#
# Usage:
#   bash human.sbatch --check
#   sbatch human.sbatch

#SBATCH --job-name=human
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

STEP_OUT=$OUT/human
RECORD_TABLE=$RECORD/human_decidual_vs_fibroblast_all_genes.csv.gz
if [ "${1:-}" = "--check" ]; then
    for f in "$HUMAN_H5AD" "$RECORD_TABLE"; do
        echo "reads:  $f ($(exists "$f"))"
    done
    echo "writes: $STEP_OUT ($(exists "$STEP_OUT"))"
    exit 0
fi

if [ -e "$STEP_OUT" ]; then
    echo "$STEP_OUT exists; skipped"
    exit 0
fi

run_notebook "$HERE/human_decidual_vs_fibroblast.ipynb" rsc "$OUT/notebooks/human_decidual_vs_fibroblast.ipynb" \
    -p OUT_DIR "$STEP_OUT" \
    -p N_CPUS "${SLURM_CPUS_PER_TASK:-8}" \
    -p H5AD "$HUMAN_H5AD" \
    -p RECORD "$RECORD_TABLE"
