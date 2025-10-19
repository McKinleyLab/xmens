#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-06:00
#SBATCH -p eddy
#SBATCH --mem=64gb
#SBATCH --array=0-1
#SBATCH -o /n/eddy_lab/users/nhilgert/decidualization/code/python/samap_integrate/10_07_logs/samap_integrate_%A_%a.out
#SBATCH -e /n/eddy_lab/users/nhilgert/decidualization/code/python/samap_integrate/10_07_logs/samap_integrate_%A_%a.err
#SBATCH --job-name=SAMap_mm_hs_pairs_0914

set -euo pipefail

echo "SLURM job started at $(date)" >&2
echo "Node: $(hostname)" >&2
echo "JOB_ID: ${SLURM_JOB_ID:-NA}  ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-NA}" >&2

# working directory
WORKDIR="/xmens/samap"
cd "$WORKDIR"

# environment
echo "Loading Python environment..." >&2
module load python
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate samap_env || { echo "Failed to activate environment 'samap_env_2'. Exiting." >&2; exit 1; }
echo "which python: $(which python)" >&2
python -V >&2

# parallelization off for now
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# inputs
# human
HUMAN_MENSTRUAL="/xmens/samap/h5ads/hs_menstrual_fibroblasts.h5ad"
HUMAN_SECRETORY="/xmens/samap/h5ads/hs_secretory_fibroblasts.h5ad"

# Mouse subsets (unchanged location)
MM_DIR="/xmens/samap/h5ads"

# Maps root (expects alignments/{hsmm|mmhs}/{hs_to_mm.txt,mm_to_hs.txt})
MAPS_DIR="$WORKDIR/alignments/"
[[ -d "$MAPS_DIR" ]] || { echo "Missing maps dir: $MAPS_DIR" >&2; exit 3; }
[[ -f "$MAPS_DIR/mmhs/hs_to_mm.txt" ]] || { echo "Missing $MAPS_DIR/mmhs/hs_to_mm.txt" >&2; exit 3; }
[[ -f "$MAPS_DIR/mmhs/mm_to_hs.txt" ]] || { echo "Missing $MAPS_DIR/mmhs/mm_to_hs.txt" >&2; exit 3; }

# --- Outputs ---
OUTDIR="$WORKDIR/samaps"
mkdir -p "$OUTDIR"
mkdir -p "$WORKDIR/logs"

# set up two pairs:
# (mm: 4_dpi,         hs: menstrual)
# (mm: 1_dpi_late,    hs: secretory)
MM_LABELS=(
  "4_dpi"         # 0
  "1_dpi_late"    # 1
)

HS_FILES=(
  "$HUMAN_MENSTRUAL"  # 0
  "$HUMAN_SECRETORY"  # 1
)

HS_TAGS=(
  "menstrual"  # 0
  "secretory"  # 1
)

# array task
IDX=${SLURM_ARRAY_TASK_ID:-0}

MM_LABEL="${MM_LABELS[$IDX]}"
DATA2="${MM_DIR}/mm_${MM_LABEL}.h5ad"   # mouse
DATA1="${HS_FILES[$IDX]}"               # human
HS_TAG="${HS_TAGS[$IDX]}"

# use the _pr.ha5d if available
export TMPDIR="${TMPDIR:-/scratch/$USER/${SLURM_JOB_ID:-samap_tmp}}"
mkdir -p "$TMPDIR" || true

FULL_D1="$WORKDIR/$DATA1"
FULL_D2="$WORKDIR/$DATA2"
[[ -f "$FULL_D1" ]] || { echo "Missing human h5ad: $FULL_D1" >&2; exit 4; }
[[ -f "$FULL_D2" ]] || { echo "Missing mouse h5ad: $FULL_D2" >&2; exit 4; }

D1_PR="${FULL_D1%.h5ad}_pr.h5ad"
D2_PR="${FULL_D2%.h5ad}_pr.h5ad"

if [[ -f "$D1_PR" ]]; then
  D1_USE="$D1_PR"
else
  D1_USE="$TMPDIR/hs_${SLURM_ARRAY_TASK_ID}.h5ad"
  cp -f "$FULL_D1" "$D1_USE"
fi

if [[ -f "$D2_PR" ]]; then
  D2_USE="$D2_PR"
else
  D2_USE="$TMPDIR/mm_${SLURM_ARRAY_TASK_ID}.h5ad"
  cp -f "$FULL_D2" "$D2_USE"
fi

echo "Run ${IDX} summary:" >&2
echo "  hs:  $D1_USE  (tag=$HS_TAG)" >&2
echo "  mm:  $D2_USE  (label=$MM_LABEL)" >&2
echo "  maps:$MAPS_DIR" >&2
echo "  out: $OUTDIR" >&2

# output filename
OUTFILE="$OUTDIR/hs_${HS_TAG}_mm_${MM_LABEL}_fibroblasts_samap.pkl"
echo "Writing to: $OUTFILE" >&2

# execute
python samap_integrate.py \
  --species1 hs \
  --species2 mm \
  --data1 "$D1_USE" \
  --data2 "$D2_USE" \
  --maps "$MAPS_DIR" \
  --out "$OUTFILE"

STATUS=$?
if [ $STATUS -eq 0 ]; then
  echo "SAMap integration (idx ${IDX}) completed successfully at $(date)." >&2
else
  echo "SAMap integration (idx ${IDX}) failed with exit code ${STATUS} at $(date)." >&2
fi

exit $STATUS