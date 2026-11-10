#!/bin/bash
# Checks the pipeline's inputs against records made when they were created (provenance/README.md):
#   labels       samples.tsv against the Samples sheet and the sequencing run written into the reads
#   FASTQs       a complete set of files per lane, and the md5 manifest once it is in sources/
#   Cell Ranger  what each lane's February run recorded: index, FASTQ folder, reference, introns, version
#   references   both references against the hashes 10x and mkref recorded when they were built
#   software     Cell Ranger, CellBender and DoubletDetection versions; sources/ against SOURCES.tsv
#   Writes: provenance/reports/check_inputs_<date>_<time>.tsv, one row per check
#   Log:    logs/check_inputs_<job ID>.out, which ends with the summary and every FAIL or PENDING row
#
# Usage, from the preprocess folder:
#   sbatch provenance/check_inputs.sbatch

#SBATCH --job-name=check_inputs
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00  # hashing the references reads about 14 GB
#SBATCH --output=logs/%x_%j.out  # job name and job ID; relative to the folder sbatch was run from

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# The pipeline folder holds config.sh and samples.tsv. Run with bash, this script sits in its provenance/ folder;
# under sbatch, SLURM runs a copy from its spool folder, and the pipeline folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")/.." && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the preprocess folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh insists on a run folder; this check writes only its report, into provenance/reports.
export OUT=$HERE/provenance/reports
source "$HERE/config.sh"

# Inputs only this check reads.
SOURCES=$HERE/provenance/sources
REF_MOUSE_INPUTS=$NH/reference/transgenic_mus_premkref  # the folder mkref built REF_MOUSE from
MOUSE_BUILD_SCRIPT=$NH/reference/mus.sh
CELLRANGER_RECORD="$NH/decidualization/revisions/velocity/{lane}/cellranger"  # the February 2026 run
QC_ENV=$(dirname "$(dirname "$PY_QC")")  # PY_QC is <environment>/bin/python

mkdir -p "$OUT"
REPORT=$OUT/check_inputs_$(date +%Y_%m_%d_%H%M).tsv

"$PY" "$HERE/provenance/check_inputs.py" \
    --samples "$HERE/samples.tsv" \
    --sheet "$SOURCES/samples_sheet_2026_11_06.csv" \
    --fastqs "$FASTQS" \
    --fastq-md5 "$SOURCES/fastq_md5_2026_11_06.txt" \
    --cellranger-record "$CELLRANGER_RECORD" \
    --ref-mouse "$REF_MOUSE" \
    --ref-mouse-inputs "$REF_MOUSE_INPUTS" \
    --transgenes "$SOURCES/transgenes/GqD_transgene.fa" "$SOURCES/transgenes/GsD_Berdeaux.fa" \
    --build-script "$MOUSE_BUILD_SCRIPT" \
    --build-script-record "$HERE/record/0_reference/mus.sh" \
    --ref-human "$REF_HUMAN" \
    --cellranger "$CELLRANGER" \
    --cellbender-env "$CELLBENDER_ENV" \
    --qc-env "$QC_ENV" \
    --sources "$SOURCES" \
    --out "$REPORT"
