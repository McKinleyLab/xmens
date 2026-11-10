#!/bin/bash
# Step 2 of 5: Cell Ranger raw matrix -> CellBender 0.2.2 remove-background, which removes ambient RNA and calls
# cells. Each array task runs one sample, task n being row n of samples.tsv.
#   Reads:  $CELLRANGER_RAW (config.sh; by default step 1's outs/raw_feature_bc_matrix.h5)
#   Writes: $OUT/cellbender/<sample>/<sample>_cellbender*; <sample>_cellbender_filtered.h5 is step 3's input
#   Logs:   $OUT/logs/2_cellbender_<sample>.log (CellBender's output), 2_cellbender_<sample>.gpu.txt (GPU, driver
#           and library versions) and logs/2_cellbender_<job ID>_<task>.out
# Command and resources are those of the July 2026 run (record/2_cellbender/01_cellbender.smk; cellbender block of
# xmens_final/workflow/profiles/slurm/config.yaml). Its jobs ran on holygpu7c0530 and holygpu7c0522, the eddy
# partition's only A40 nodes, so no nodes are excluded. A sample whose filtered .h5 exists is skipped.
#
# Usage:
#   bash 2_cellbender.sbatch --check        print each sample's input and command; run nothing
#   sbatch 2_cellbender.sbatch              all samples
#   sbatch --array=3 2_cellbender.sbatch    row 3 only

#SBATCH --job-name=2_cellbender
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --gres=gpu:nvidia_a40:1  # one A40 GPU
#SBATCH --cpus-per-task=1  # the work runs on the GPU
#SBATCH --mem=64G
#SBATCH --time=08:00:00  # samples took 1.2-1.8 h in July 2026
#SBATCH --array=1-15  # one task per row of samples.tsv
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

# config.sh sets $OUT, $CELLRANGER_RAW, $CONDA, $CELLBENDER_ENV, N_SAMPLES, and field <row> <column>, which prints
# one value from samples.tsv.
source "$HERE/config.sh"

# Sets SAMPLE, RAW, DIR and COMMAND for row $1 of samples.tsv. RAW is CELLRANGER_RAW with {lane} replaced by the
# row's lane (bash ${VAR//pattern/replacement}). COMMAND is a bash array, so it can be printed and run without
# paths being split at spaces.
sample_command() {
    SAMPLE=$(field "$1" sample_id)
    RAW=${CELLRANGER_RAW//\{lane\}/$(field "$1" lane)}
    DIR=$OUT/cellbender/$SAMPLE

    # --total-droplets-included: how many barcodes, ranked by UMI count, CellBender models; the rest are taken as
    # empty droplets. --cuda runs it on the GPU. Everything else is CellBender 0.2.2's default (150 epochs,
    # false-positive rate 0.01).
    # The environment's cellbender launcher names its interpreter by the old /n/eddy_lab mount, which no longer
    # exists, so the launcher is run by the environment's own python: same program, same interpreter.
    COMMAND=(
        "$CELLBENDER_ENV/bin/python" "$CELLBENDER_ENV/bin/cellbender" remove-background
        --input "$RAW"
        --output "$DIR/${SAMPLE}_cellbender.h5"
        --total-droplets-included 50000
        --cuda
    )
}

# --check: for every row, print whether the input exists, whether the sample is done, and the command; then stop
# without running anything.
if [ "${1:-}" = "--check" ]; then
    for row in $(seq 1 "$N_SAMPLES"); do
        sample_command "$row"

        input="input missing"
        if [ -e "$RAW" ]; then
            input="input found"
        fi

        state="to run"
        if [ -e "$DIR/${SAMPLE}_cellbender_filtered.h5" ]; then
            state="done"
        fi

        echo "$row $SAMPLE: $input, $state"
        echo "    ${COMMAND[*]}"
    done
    exit 0
fi

# The array task ID, which sbatch sets, is the row of samples.tsv this task runs; without one, the script stops
# with the message below.
ROW=${SLURM_ARRAY_TASK_ID:?submit with sbatch, or set SLURM_ARRAY_TASK_ID to a row of samples.tsv}
sample_command "$ROW"

if [ -e "$DIR/${SAMPLE}_cellbender_filtered.h5" ]; then
    echo "$SAMPLE: $DIR/${SAMPLE}_cellbender_filtered.h5 exists; nothing to do"
    exit 0
fi
if [ ! -e "$RAW" ]; then
    echo "$SAMPLE: input $RAW not found" >&2
    exit 1
fi
mkdir -p "$DIR" "$OUT/logs"

# A batch shell does not inherit conda, so conda goes on PATH before the environment is activated.
export PATH=$CONDA/bin:$PATH
source "$CONDA/etc/profile.d/conda.sh"
conda activate "$CELLBENDER_ENV"

# The GPU and library versions this sample ran on, next to its log. CellBender trains on the GPU; two runs on the
# same node and versions gave identical output, so these are what a difference would be traced to.
# Informational only, so a failure here never stops the run.
{
    echo "host $(hostname), CUDA_VISIBLE_DEVICES ${CUDA_VISIBLE_DEVICES:-unset}"
    nvidia-smi --query-gpu=name,uuid,driver_version,memory.total --format=csv,noheader || true
    ls -d "$CELLBENDER_ENV"/lib/python3*/site-packages/cellbender-*.dist-info 2>/dev/null | sed 's|.*/||' || true
    python -c "import torch; print('torch', torch.__version__, '| CUDA', torch.version.cuda,
        '| cuDNN', torch.backends.cudnn.version(), '| GPU', torch.cuda.get_device_name(0))" || true
} > "$OUT/logs/2_cellbender_$SAMPLE.gpu.txt" 2>&1

# CellBender's progress goes to the step log; the SLURM log gets the start and end times.
echo "$SAMPLE: started $(date) on $(hostname)"
"${COMMAND[@]}" > "$OUT/logs/2_cellbender_$SAMPLE.log" 2>&1
echo "$SAMPLE: finished $(date)"
