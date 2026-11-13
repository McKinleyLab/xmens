bash <<'SCRIPT'
# Reproducibility check of the mouse pseudotime of record (10 Nov): the rerun's code and object against the final
# ones; appends the result to RUNS.md.
set -euo pipefail
NH=/n/lab_storage/eddy_lab/users/nhilgert
NS=/n/netscratch/eddy_lab/Everyone/nhilgert
cd $NH/xmens
FINAL=/n/lab_storage/eddy_lab/Lab/mckinley/data/20261110_mouse_only_eFib_dec_final_nick_pdt.h5ad
FINAL_SHA=e31372ba5fdd706ad620e8259e34216bc96b8f3b6394e4d2eedc0baeec8760dd   # pinned in mi/ and samap/config.sh
RERUN=/n/lab_storage/eddy_lab/Lab/mckinley/data/20261110_mouse_only_eFib_dec_final_nick_pdt_rerun.h5ad
CHECK=$NS/2026_11_12_mouse_pdt_rerun_check                                   # only holds the h5diff log
mkdir -p $CHECK
H5DIFF=$NH/.conda/envs/cellbender_022/bin/h5diff
PY=$NH/.uv/envs/rsc_uv/bin/python
NB_FINAL=pseudotime/pseudotime_mouse_final_1110.ipynb
NB_RERUN_SRC=pseudotime/pseudotime_mouse_final_1110_rerun.ipynb
NB_RERUN=pseudotime/record/pseudotime_mouse_final_1110_rerun.executed_20261112.ipynb
test ! -e $NB_RERUN && cp -p $NB_RERUN_SRC $NB_RERUN

# 0. the final object is still the one the configs pin
final_sha=$(sha256sum "$FINAL" | cut -c1-64)
[ "$final_sha" = "$FINAL_SHA" ] && final_note="sha256 as pinned" || final_note="sha256 $final_sha, NOT the pinned value"

# 1. code cells of the rerun vs the committed notebook (expect only the parameters cell, [0])
code_diff=$($PY - "$NB_FINAL" "$NB_RERUN" <<'EOF'
import json, sys
code = lambda p: ["".join(c["source"]) for c in json.load(open(p))["cells"] if c["cell_type"] == "code"]
a, b = code(sys.argv[1]), code(sys.argv[2])
print(f"{len(a)} vs {len(b)} code cells; differing code cells {[i for i, (x, y) in enumerate(zip(a, b)) if x != y]}")
EOF
)

# 2. the two objects, everything except the write time; on a difference, list where
set +e
$H5DIFF --exclude-path /uns/pseudotime/computed "$FINAL" "$RERUN" > $CHECK/h5diff.txt 2>&1
code=$?
set -e
if [ $code -eq 0 ]; then
    verdict="identical (h5diff exit 0)"
else
    paths=$(grep -oE '^(dataset|attribute|group): <[^>]+>' $CHECK/h5diff.txt | sort -u | head -20 | tr '\n' ' ' || true)
    verdict="DIFFERENT (h5diff exit $code): ${paths:-see $CHECK/h5diff.txt}"
fi
written=$($PY - "$RERUN" <<'EOF'
import h5py, sys
g = h5py.File(sys.argv[1], "r")["uns/pseudotime"]
v = g["computed"][()] if "computed" in g else "n/a"
print(v.decode() if isinstance(v, bytes) else v)
EOF
)
sha=$(sha256sum "$RERUN" | cut -c1-64)

cat >> pseudotime/record/RUNS.md <<EOT

## $(date +%F) — mouse pseudotime, rerun of the 10 Nov object of record (reproducibility check)
Code:     $NB_RERUN_SRC, a copy of $NB_FINAL ($code_diff)
Output:   $RERUN  sha256 $sha  (written $written)
Check:    h5diff vs $FINAL ($final_note), excluding /uns/pseudotime/computed: $verdict
Executed: $NB_RERUN
Status:   check
EOT
tail -n 7 pseudotime/record/RUNS.md
SCRIPT