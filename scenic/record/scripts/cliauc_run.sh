#!/bin/bash
# Re-score the floor-6 consensus regulons with pySCENIC's own command-line AUCell, unmodified, as in the protocol:
#   pyscenic aucell <filtered loom> <regulons.yaml> --output <csv> [--weights] --seed 0
# Input loom = {sp}_dsc_pooled.loom (all cells, raw counts, genes detected in >= 1 % of cells). pySCENIC 0.12.1's CLI
# cannot read pickled (.dat) regulons (ctxcore openfile() passes an encoding to a binary open), so the regulons are
# written to YAML with pySCENIC's own save_to_yaml and verified to round-trip exactly with load_from_yaml.
# The env moved from /n/eddy_lab to /n/lab_storage; bin/pyscenic keeps the old shebang, so the same entry script is run
# with the env's own interpreter. Writes only new files; refuses if the output folder exists.
set -euo pipefail
set -o noclobber
LANE=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/scenic_v2_2026_08_28
ENV=/n/lab_storage/eddy_lab/users/nhilgert/.conda/envs/pyscenic
PYR=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python
NEW=consensus_sketch15k_weighted_soft_occ6_pyscenic_cli
FAILED=consensus_sketch15k_weighted_soft_occ6_cliaucell
export TMPDIR=/n/lab_storage/eddy_lab/users/nhilgert/tmp XDG_CACHE_HOME=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch
export NUMBA_CACHE_DIR=$XDG_CACHE_HOME/numba MPLCONFIGDIR=$XDG_CACHE_HOME/mpl
for SP in mouse human; do
  if [ -e "$LANE/scenic_seeds/$SP/$NEW" ]; then echo "EXISTS: $SP $NEW -- refusing to overwrite"; exit 3; fi
done
# note in the folder left by the failed 2026-10-26 attempt (nothing removed)
if [ -d "$LANE/scenic_seeds/mouse/$FAILED" ] && [ ! -e "$LANE/scenic_seeds/mouse/$FAILED/NOTE.md" ]; then
  cat > "$LANE/scenic_seeds/mouse/$FAILED/NOTE.md" <<EOF
# Failed attempt, 2026-10-26 11:09 -- not a result
pySCENIC 0.12.1's command line could not read the pickled regulons (consensus_regulons.dat, a byte-identical copy of
consensus_sketch15k_weighted_soft_occ6/consensus_regulons.pkl) and stopped; auc_mtx_unweighted.csv is the empty file it
had opened. The completed run is in ../$NEW. Kept, not deleted; safe to remove.
EOF
fi
echo "start $(date '+%F %T')" > run.log
for SP in mouse human; do
  SRC=$LANE/scenic_seeds/$SP/consensus_sketch15k_weighted_soft_occ6
  OUT=$LANE/scenic_seeds/$SP/$NEW
  mkdir "$OUT"
  "$ENV/bin/python" - "$SRC/consensus_regulons.pkl" "$OUT/consensus_regulons.yaml" <<'EOF' | tee -a run.log
import pickle, sys
from pyscenic.utils import save_to_yaml, load_from_yaml
src, dst = sys.argv[1], sys.argv[2]
regs = pickle.load(open(src, "rb"))
save_to_yaml(regs, dst)
back = load_from_yaml(dst)
assert len(back) == len(regs)
a = {r.name: r for r in regs}; b = {r.name: r for r in back}
assert set(a) == set(b)
for n in a:
    assert dict(a[n].gene2weight) == dict(b[n].gene2weight), n
    assert a[n].transcription_factor == b[n].transcription_factor, n
print(f"yaml round-trip exact: {len(regs)} regulons, {sum(len(r.genes) for r in regs)} TF-target pairs, weights identical")
EOF
  LOOM=$LANE/scenic_seeds/$SP/${SP}_dsc_pooled.loom
  for W in unweighted weighted; do
    FLAG=""; [ "$W" = weighted ] && FLAG="--weights"
    CMD="$ENV/bin/python $ENV/bin/pyscenic aucell $LOOM $OUT/consensus_regulons.yaml --output $OUT/auc_mtx_${W}.csv --num_workers 8 --seed 0 $FLAG"
    echo "[$(date '+%F %T')] $SP $W: $CMD" | tee -a run.log
    $CMD > "aucell_${SP}_${W}.log" 2>&1
    test -s "$OUT/auc_mtx_${W}.csv"
    echo "[$(date '+%F %T')] $SP $W done" | tee -a run.log
  done
  {
    echo "# $NEW"
    echo
    echo "Floor-6 consensus regulons of consensus_sketch15k_weighted_soft_occ6 (same TFs, targets and occurrence weights),"
    echo "scored with pySCENIC 0.12.1's command-line AUCell, unmodified, as in the published protocol, on"
    echo "${SP}_dsc_pooled.loom: all cells, raw counts, genes detected in >= 1 % of cells (the network-inference gene set)."
    echo "Made $(date '+%F %T') by cliauc_run.sh. Commands:"
    grep " $SP .*: " run.log | sed 's/^/    /'
    echo
    echo "- consensus_regulons.yaml: the regulons, written with pyscenic.utils.save_to_yaml; verified to reload exactly."
    echo "- auc_mtx_weighted.csv: PRIMARY. --weights, i.e. each target weighted by its occurrence fraction across runs."
    echo "- auc_mtx_unweighted.csv: CLI defaults (no --weights), for comparison."
    echo "Both: auc_threshold 0.05 (default), ranking seed 0 (the protocol leaves it random)."
  } > "$OUT/README.md"
  sha256sum "$OUT"/* > "sha_$SP.txt"
done
"$PYR" cliauc_check.py 2>&1 | tee check.log
