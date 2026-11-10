#!/bin/bash
#SBATCH -p eddy
#SBATCH -t 0-00:40
#SBATCH --mem=96G
#SBATCH -c 4
#SBATCH -J prep_raw
#SBATCH -o /n/lab_storage/eddy_lab/users/nhilgert/code/samap-scvi/full_fibdec/logs/%x_%j.out
# Build RAW-COUNT SAMap inputs (same gene sets/obs as the CP10K inputs; only .X differs).
set -uo pipefail
DIR=/n/lab_storage/eddy_lab/users/nhilgert/code/samap-scvi/full_fibdec
PY=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/samap_harmony_uv/bin/python
NS=/n/netscratch/eddy_lab/Everyone/nhilgert
module load gcc/14.2.0-fasrc01
cd "$DIR"
mkdir -p inputs_raw logs

"$PY" prep_full_input.py --scale raw \
    --src "$NS/2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad" \
    --species hs --out inputs_raw/hs_fib_dec_all_raw.h5ad || exit 1
"$PY" prep_full_input.py --scale raw \
    --src "$NS/2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad" \
    --species mm --out inputs_raw/mm_fib_dec_all_raw.h5ad || exit 1
echo "=== prep done ($(date)) ==="
