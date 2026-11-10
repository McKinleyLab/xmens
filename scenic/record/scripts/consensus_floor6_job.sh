#!/bin/bash
# Floor-6 consensus variant of the SCENIC v2 lane (targets recurring in >= 6 runs, i.e. "more than 5" as in Ni et al. and
# the VSN multi-run code behind the Fly Cell Atlas). Stage 3 only: consensus + AUCell from the existing 100 per-run regulon
# pickles; no GRNBoost2 / cisTarget. Same flags as the floor-5 and floor-80 variants except --min-occurrence.
# Writes new directories consensus_sketch15k_weighted_soft_occ6/ in the lane; refuses to overwrite.
set -euo pipefail
REV=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/scenic_v2_2026_08_28
PY=/n/lab_storage/eddy_lab/users/nhilgert/.conda/envs/pyscenic/bin/python
ENV=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv
FIG=/n/netscratch/eddy_lab/Everyone/nhilgert/figures
C=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch
export XDG_CACHE_HOME=$C NUMBA_CACHE_DIR=$C/numba MPLCONFIGDIR=$C/mpl TMPDIR=/n/lab_storage/eddy_lab/users/nhilgert/tmp
export IPYTHONDIR=$C/ipython JUPYTER_RUNTIME_DIR=$C/jupyter_runtime JUPYTER_DATA_DIR=$C/jupyter_data
mkdir -p "$IPYTHONDIR" "$JUPYTER_RUNTIME_DIR" "$JUPYTER_DATA_DIR"
for SP in mouse human; do
  if [ -e "$REV/scenic_seeds/$SP/consensus_sketch15k_weighted_soft_occ6" ]; then echo "EXISTS: $SP occ6 -- refusing to overwrite"; exit 3; fi
done

# 0. smoke test: the netscratch notebook runs with its csv next to it (run on a copy; netscratch is not written)
set +e
mkdir -p nbsmoke
cp "$FIG/scenic/fig_regulon_quadrants_floor5.ipynb" "$FIG/scenic/quadfig.py" "$FIG/scenic/regulon_quadrants_floor5_input.csv" nbsmoke/
"$ENV/bin/jupyter" nbconvert --to notebook --execute nbsmoke/fig_regulon_quadrants_floor5.ipynb --output smoke_executed.ipynb > nbsmoke.log 2>&1
echo "nbconvert exit $?" > nbsmoke.txt
grep -o '111 TFs from [^;]*' nbsmoke/smoke_executed.ipynb | head -1 >> nbsmoke.txt
ls nbsmoke >> nbsmoke.txt
set -e

# 1. floor-6 consensus + AUCell, both species
for SP in mouse human; do
  ROOT=$REV/scenic_seeds/$SP
  OUT=$ROOT/consensus_sketch15k_weighted_soft_occ6
  "$PY" "$REV/scripts/build_consensus_aucell.py" --species $SP --rev "$REV" --h5ad "$ROOT/${SP}_dsc_pooled.h5ad" \
      --timepoint-col timepoint --n-cpu 8 --sketch 15k --mode weighted_soft --min-occurrence 6 \
      --out-dir "$OUT" > consensus_$SP.log 2>&1
  cp "$OUT/recurrence_summary.txt" summary_occ6_$SP.txt
  cp "$OUT/recurrence_report.tsv" report_occ6_$SP.tsv
done

# 2. rho of regulon activity with pseudotime, floors 5 and 6 on the same obs
"$PY" auc_vs_pt_floor6.py > auc_vs_pt.log 2>&1

# 3. floor 6 vs floor 5
"$PY" compare_occ6.py > compare.log 2>&1
cat compare.log
