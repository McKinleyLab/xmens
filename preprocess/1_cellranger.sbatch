#!/bin/bash
# Step 1 of 5: FASTQs -> Cell Ranger 7.2.0 count matrices. Each array task runs one lane, task n being row n of
# samples.tsv.
#   Reads:  the lane's FASTQs and the species' reference (paths in config.sh)
#   Writes: $OUT/cellranger/<lane>/; its outs/raw_feature_bc_matrix.h5 is step 2's input
#   Log:    logs/1_cellranger_<job ID>_<task>.out
# A lane whose raw matrix exists is skipped. Rerunning a failed lane resumes it; if Cell Ranger reports that the
# pipestance is locked, the killed run left $OUT/cellranger/<lane>/_lock behind, and that file has to go first.
#
# Usage:
#   bash 1_cellranger.sbatch --check                  print each lane's command and FASTQ count; run nothing
#   sbatch 1_cellranger.sbatch                        all lanes
#   sbatch --array=3 1_cellranger.sbatch              row 3 only
#   SLURM_ARRAY_TASK_ID=3 bash 1_cellranger.sbatch    row 3 in the current shell (an interactive session)

#SBATCH --job-name=1_cellranger
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=8  # Cell Ranger is told to use 8 cores (--localcores below)
#SBATCH --mem=100G  # Cell Ranger is held to 96 GB (--localmem below), which leaves room for the rest of the job
#SBATCH --time=2-00:00:00  # 2 days; lanes took 2 to 6+ hours in November 2026
#SBATCH --array=1-15  # one task per row of samples.tsv
#SBATCH --exclude=holygpu8a26504,holygpu8a26505,holygpu8a26506
#SBATCH --output=logs/%x_%A_%a.out  # job name, job ID, task; relative to the folder sbatch was run from

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

# config.sh sets the paths ($OUT, $CELLRANGER, $FASTQS, $REF_MOUSE, $REF_HUMAN, $CACHE), N_SAMPLES, and
# field <row> <column>, which prints one value from samples.tsv.
source "$HERE/config.sh"

# Sets LANE, FASTQ_DIR, PREFIX, REF and COMMAND for row $1 of samples.tsv. COMMAND is a bash array, one word per
# element, so it can be printed (--check) and run ("${COMMAND[@]}") without paths being split at spaces.
lane_command() {
    LANE=$(field "$1" lane)
    FASTQ_DIR=$FASTQS/$(field "$1" fastq_dir)
    PREFIX=$(field "$1" fastq_prefix)

    REF=$REF_MOUSE
    if [ "$(field "$1" species)" = human ]; then
        REF=$REF_HUMAN
    fi

    # The February 2026 command (record/1_cellranger/feb_cmdlines/) plus --no-bam, which skips the BAM and leaves
    # the matrices unchanged (byte-identical to February in the November 2026 re-run), and --disable-ui, which
    # turns off Cell Ranger's web page. --sample is the FASTQ file-name prefix (SI-TT-A5 for
    # SI-TT-A5_S1_L001_R1_001.fastq.gz); --localcores and --localmem (GB) keep Cell Ranger within the SLURM request.
    COMMAND=(
        "$CELLRANGER" count
        --id="$LANE"
        --transcriptome="$REF"
        --fastqs="$FASTQ_DIR"
        --sample="$PREFIX"
        --output-dir="$OUT/cellranger/$LANE"
        --localcores=8
        --localmem=96
        --no-bam
        --disable-ui
    )
}

# --check: for every row, print the number of R1 FASTQ files found, the reference, whether the lane is done, and
# the command; then stop without running anything.
if [ "${1:-}" = "--check" ]; then
    for row in $(seq 1 "$N_SAMPLES"); do
        lane_command "$row"

        # 0 files means a wrong FASTQ folder or prefix. ls fails when nothing matches, so 2>/dev/null and
        # || true let --check report 0 instead of stopping.
        n_r1=$(ls "$FASTQ_DIR/$PREFIX"_S*_R1_001.fastq.gz 2>/dev/null | wc -l | tr -d ' ' || true)

        state="to run"
        if [ -e "$OUT/cellranger/$LANE/outs/raw_feature_bc_matrix.h5" ]; then
            state="done"
        fi

        # ${REF##*/} is the reference's folder name; ${COMMAND[*]} joins the command into one line.
        echo "$row $LANE: $n_r1 R1 FASTQ files, reference ${REF##*/}, $state"
        echo "    ${COMMAND[*]}"
    done
    exit 0
fi

# The array task ID, which sbatch sets, is the row of samples.tsv this task runs; without one, the script stops
# with the message below.
ROW=${SLURM_ARRAY_TASK_ID:?submit with sbatch, or set SLURM_ARRAY_TASK_ID to a row of samples.tsv}
lane_command "$ROW"

if [ -e "$OUT/cellranger/$LANE/outs/raw_feature_bc_matrix.h5" ]; then
    echo "$LANE: $OUT/cellranger/$LANE/outs/raw_feature_bc_matrix.h5 exists; nothing to do"
    exit 0
fi

# Cell Ranger writes temporary files and a __<lane>.mro into its working folder, so both stay in the run folder.
# Its cache and settings would otherwise go to the home directory.
mkdir -p "$OUT/cellranger/tmp/$LANE" "$CACHE/cellranger"
export TMPDIR=$OUT/cellranger/tmp/$LANE
export XDG_CACHE_HOME=$CACHE/cellranger
export XDG_CONFIG_HOME=$CACHE/cellranger
cd "$OUT/cellranger"

echo "$LANE: started $(date) on $(hostname)"
"${COMMAND[@]}"
echo "$LANE: finished $(date)"
