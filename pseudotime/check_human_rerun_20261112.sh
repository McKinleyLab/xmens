#!/bin/bash
# Reproducibility check of the final human pseudotime: the rerun's code and object against the final ones;
# appends the result to RUNS.md.
set -euo pipefail
NH=/n/lab_storage/eddy_lab/users/nhilgert
cd $NH/xmens
FINAL=/n/lab_storage/eddy_lab/Lab/mckinley/data/20261112_Human_fibdec_final_nick_pdt.h5ad
RERUN=$(ls /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_human_reroot_rerun/*.h5ad)
H5DIFF=$NH/.conda/envs/cellbender_022/bin/h5diff
PY=$NH/.uv/envs/rsc_uv/bin/python
NB_FINAL=pseudotime/pseudotime_human_final_1112.ipynb
NB_RERUN_SRC=pseudotime/pseudotime_human_final_1112_rerun.ipynb
NB_RERUN=pseudotime/record/pseudotime_human_final_1112_rerun.executed_20261112.ipynb
test ! -e $NB_RERUN && cp -p $NB_RERUN_SRC $NB_RERUN

# 1. code cells of the rerun vs the final notebook (expect only the parameters cell, 0, to differ)
code_diff=$($PY - "$NB_FINAL" "$NB_RERUN" <<'EOF'
import json, sys
code = lambda p: ["".join(c["source"]) for c in json.load(open(p))["cells"] if c["cell_type"] == "code"]
a, b = code(sys.argv[1]), code(sys.argv[2])
print(f"{len(a)} vs {len(b)} code cells; differing code cells {[i for i, (x, y) in enumerate(zip(a, b)) if x != y]}")
EOF
)
# 2. the two objects, everything except uns['pseudotime']['computed'] (the write time)
if $H5DIFF --exclude-path /uns/pseudotime/computed "$FINAL" "$RERUN" > /dev/null; then
    verdict="identical (h5diff exit 0)"
else
    verdict="DIFFERENT (h5diff exit $?); rerun h5diff without > /dev/null to see where"
fi
written=$($PY -c "import h5py, sys; v = h5py.File(sys.argv[1], 'r')['uns/pseudotime/computed'][()]; \
print(v.decode() if isinstance(v, bytes) else v)" "$RERUN")
sha=$(sha256sum "$RERUN" | cut -c1-64)

cat >> pseudotime/record/RUNS.md <<EOT

## $(date +%F) — human pseudotime, rerun of the final (reproducibility check)
Code:     $NB_RERUN_SRC, a copy of $NB_FINAL ($code_diff)
Output:   $RERUN  sha256 $sha  (written $written)
Check:    h5diff vs $FINAL, excluding /uns/pseudotime/computed: $verdict
Executed: $NB_RERUN
Status:   check
EOT
tail -n 8 pseudotime/record/RUNS.md