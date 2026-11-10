#!/bin/bash
#SBATCH -p eddy
#SBATCH -t 0-10:00
#SBATCH --mem=128G
#SBATCH -c 8
#SBATCH -J samap_raw
#SBATCH -o /n/lab_storage/eddy_lab/users/nhilgert/code/samap-scvi/full_fibdec/logs/%x_%j.out
#
# SAMap on RAW COUNTS — the input the paper specifies. Same objects, same genes,
# same cells, same SAMAP settings as the CP10K experiment; only .X scale differs.
#
#   vanilla_raw   : stock SAM, min_expression=1  (full by-the-book path)
#   harmony_raw   : + harmonypy inside each SAM iteration (batch_key='sample')
#   vanilla_raw_nofilter : min_expression=0, isolates the low-count filter from
#                          the normalization target (cell_median vs CP10K's 1e4)
set -uo pipefail

DIR=/n/lab_storage/eddy_lab/users/nhilgert/code/samap-scvi/full_fibdec
PY=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/samap_harmony_uv/bin/python
MAPS=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/samap_integrate/homology/hsmm_hmmer
HS=inputs_raw/hs_fib_dec_all_raw.h5ad
MM=inputs_raw/mm_fib_dec_all_raw.h5ad

module load gcc/14.2.0-fasrc01
cd "$DIR"
export PYTHONPATH="$DIR"
mkdir -p runs_raw logs

for f in "$HS" "$MM"; do [ -f "$f" ] || { echo "missing $f"; exit 1; }; done

run () {  # $1=out name  $2=arm  $3=min_expression
    echo "=== arm: $1 (arm=$2 min_expression=$3) ($(date)) ==="
    "$PY" run_arm.py \
        --hs "$HS" --mm "$MM" --maps-dir "$MAPS" \
        --arm "$2" --out-dir "runs_raw/$1" \
        --batch-key sample --min-expression "$3" \
        --n-iter 3 --ncpus "$SLURM_CPUS_PER_TASK" \
        || echo "!!! ARM $1 FAILED (continuing) !!!"
}

run vanilla_raw          vanilla 1
run harmony_raw          harmony 1
run vanilla_raw_nofilter vanilla 0

echo "=== compare ($(date)) ==="
"$PY" compare_arms.py runs_raw/vanilla_raw runs_raw/harmony_raw \
    runs_raw/vanilla_raw_nofilter | tee runs_raw/summary.txt
echo "=== done ($(date)) ==="
