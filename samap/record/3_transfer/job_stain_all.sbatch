#!/bin/bash
#SBATCH -p eddy
#SBATCH -t 0-03:00
#SBATCH --mem=64G
#SBATCH -c 4
#SBATCH -J stain_all
#SBATCH -o /n/lab_storage/eddy_lab/users/nhilgert/code/samap-scvi/full_fibdec/logs/%x_%j.out
# Cross-species label-transfer stainings, coarse AND fine, both directions,
# for the CP10K arms and the raw-count arms.
set -uo pipefail
DIR=/n/lab_storage/eddy_lab/users/nhilgert/code/samap-scvi/full_fibdec
PY=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/samap_harmony_uv/bin/python
module load gcc/14.2.0-fasrc01
cd "$DIR"

for RUN in runs/vanilla runs/harmony \
           runs_raw/vanilla_raw runs_raw/harmony_raw runs_raw/vanilla_raw_nofilter; do
  if [ -f "$RUN/samap.pkl" ]; then
    echo "=== stain $RUN ($(date)) ==="
    "$PY" label_stain_full.py --samap-pkl "$RUN/samap.pkl" \
        --out-dir "$RUN/staining" || echo "!!! $RUN FAILED (continuing) !!!"
  else
    echo "=== skip $RUN (no pkl) ==="
  fi
done
echo "=== done ($(date)) ==="
