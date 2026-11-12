# The CellRank row repair in the pseudotime: where it came from, and that no object of record used it

Checked 12 Nov 2026. Each claim is followed by the command that shows it and its output as printed that day.

**What the problem is.** CellRank's RealTimeKernel builds each nucleus's row of the transition matrix from 0.2 x its
same-day neighbors and 0.8 x its moscot transport row to the next timepoint, and refuses the whole matrix if any row
does not sum to 1 (`ValueError: Transition matrix is not row stochastic`). The transport row comes from entropic optimal
transport (epsilon = 0.001), whose entries scale like exp(-distance / epsilon). For a nucleus far from every nucleus of
the next timepoint, the row is either exactly 0, which moscot replaces with a uniform row ("empty row(s) ... Using
uniform distribution", in every run), or tiny but nonzero, which turns into infinite entries when normalized. CellRank
rejects the second kind, and nothing upstream handles it.

**The repair** replaced each such row (uniform over the nucleus's own timepoint, or, in the November version, its
same-day transitions only), renormalized every row, and handed the matrix back to CellRank, which checked it again.

**Where it ran.** It fired only in exploratory runs: about 20 Apr 2026 (4 rows), 28 Aug 2026 (5 rows; then 4 in the
patched version), and one interactive run on 7 Nov 2026 (5 rows). It never fired in a run that wrote an object of
record: the human run of record on the 2025 objects and the Nov 6 re-runs logged 0, and the three Nov 10 runs that wrote
the `_nick_pdt` objects ran without the repair, with CellRank accepting every matrix. It was removed from the November
notebooks; a run that hits the case now stops.

**Correction to what we said on 10-11 Nov.** "It never fired in any November run" is too broad: the interactive run of
7 Nov repaired 5 rows (step 10). The accurate statement is the one in the notebooks' comment: no run of record hit it.
For the mouse run of record on the 2025 objects (22 Oct), no repaired-rows count was found; the Nov 6 re-run on the same
object logged 0.

## Setup

```bash
NH=/n/lab_storage/eddy_lab/users/nhilgert
R=$NH/decidualization/revisions/cellrank
D=/n/lab_storage/eddy_lab/Lab/mckinley/data
cell_out() {  # saved output of cell $2 of notebook $1 ("-" reads the notebook from stdin)
  python3 -c 'import json, sys
nb = json.load(sys.stdin if sys.argv[1] == "-" else open(sys.argv[1]))
for o in nb["cells"][int(sys.argv[2])].get("outputs", []):
    if "text" in o: print("".join(o["text"]), end="")
    elif o.get("output_type") == "error": print(o["ename"] + ": " + o["evalue"])' "$1" "$2"
}
cell_src() {  # first $3 lines of the source of cell $2 of notebook $1 ("-" reads stdin)
  python3 -c 'import json, sys
nb = json.load(sys.stdin if sys.argv[1] == "-" else open(sys.argv[1]))
print("\n".join("".join(nb["cells"][int(sys.argv[2])]["source"]).splitlines()[:int(sys.argv[3])]))' "$1" "$2" "$3"
}
```

## 1

The first failure, about 20 Apr 2026: an exploratory mouse run (33,717 nuclei, no uninduced) in `decidualization/revisions/cellrank/realtime_informed_pseudotime_variants.ipynb`. moscot replaces its exactly empty rows with a uniform row; CellRank still rejects 4 rows.

```bash
ls -l --time-style=+%F_%H:%M $R/realtime_informed_pseudotime_variants.ipynb | awk '{print $6, $7}'
cell_out $R/realtime_informed_pseudotime_variants.ipynb 25
```
```
2026-04-20_16:10 /n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/cellrank/realtime_informed_pseudotime_variants.ipynb
Using automatic `threshold=0.0`
WARNING: Coupling at `(1.0, 2.0)` contains `4` empty row(s), e.g., due to unbalancedness or 0s in the source marginals. Using uniform distribution
WARNING: Coupling at `(3.0, 4.0)` contains `4` empty row(s), e.g., due to unbalancedness or 0s in the source marginals. Using uniform distribution
ValueError: Transition matrix is not row stochastic, 4 rows do not sum to 1.
```

## 2

The diagnosis (same notebook, next cell): the matrix rebuilt from the couplings has rows that sum to infinity. These are transport rows that underflowed to a tiny, nonzero sum and blew up when normalized. (The 24,536 rows "non-stochastic" at rtol 1e-12 are float32 rounding at the 1e-8 level; CellRank's 4 are the infinite ones.)

```bash
cell_out $R/realtime_informed_pseudotime_variants.ipynb 26
```
```
Shape: (33717, 33717)
Row sums — min: 1.000000
Row sums — max: inf
Row sums — mean: inf
Row sums — median: 1.000000

Non-stochastic rows: 24536
Bad row sums — min: 1.000000
Bad row sums — max: inf
Bad row sum distribution: [0.9999998452 0.9999999881 1.0000000119 1.0000000315          nan]
```

## 3

The first repair (next cell): each row with a zero or non-finite sum becomes uniform over the nuclei of its own timepoint; every row is renormalized; the matrix is handed back to CellRank, which checks it again.

```bash
cell_out $R/realtime_informed_pseudotime_variants.ipynb 27
```
```
Replacing 4 non-finite/zero rows with within-timepoint uniform...
After fix — row sums: min=0.999999999999999, max=1.000000000000001
Non-stochastic rows (rtol=1e-12): 0
Non-stochastic rows (rtol=1e-6):  0

rtk.transition_matrix set successfully.
```

## 4

The same cells copied into the human notebook (28 May 2026, cells 27-30): CellRank accepts the human matrix directly and the repair finds nothing to do.

```bash
ls -l --time-style=+%F_%H:%M $R/realtime_informed_pseudotime_human.ipynb | awk '{print $6, $7}'
cell_out $R/realtime_informed_pseudotime_human.ipynb 27
cell_out $R/realtime_informed_pseudotime_human.ipynb 29 | head -1
```
```
2026-05-28_11:43 /n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/cellrank/realtime_informed_pseudotime_human.ipynb
Using automatic `threshold=2.147570896515066e-19`
Transition matrix: csr_matrix  shape=(44756, 44756)
Row sums:  min=1.000000  max=1.000000
Replacing 0 non-finite/zero rows with within-timepoint uniform...
```

## 5

The second failure, 28 Aug 2026 (time stamp of the run's GPU log line): mouse on the 2026_08_28 objects, all 43,197 fib/dec nuclei including the apoptotic ones, with the same moscot and CellRank settings as every later run. 5 rows rejected. `xmens_final/notebooks/fibdec_pseudotime.ipynb`, in git from checkpoint `b9fa19c`.

```bash
git -C $NH/xmens_final show b9fa19c:notebooks/fibdec_pseudotime.ipynb | cell_out - 13 | head -1
git -C $NH/xmens_final show b9fa19c:notebooks/fibdec_pseudotime.ipynb | cell_out - 6 | grep -E 'variant_name|remove_apoptotic|moscot_|rtk_'
git -C $NH/xmens_final show b9fa19c:notebooks/fibdec_pseudotime.ipynb | cell_out - 28 | grep -m1 -o 'E0828 [0-9:.]*'
git -C $NH/xmens_final show b9fa19c:notebooks/fibdec_pseudotime.ipynb | cell_out - 30
```
```
AnnData object with n_obs × n_vars = 43197 × 33698
  variant_name                   = v1_all_cells
  remove_apoptotic               = False
  moscot_epsilon                 = 0.001
  moscot_tau_a                   = 0.95
  moscot_tau_b                   = 0.95
  rtk_conn_weight                = 0.2
  rtk_self_transitions           = all
  rtk_threshold                  = auto
  moscot_gene_set                = mouse
E0828 22:44:55.419984
Using automatic `threshold=0.0`
WARNING: Coupling at `(3.0, 4.0)` contains `1` empty row(s), e.g., due to unbalancedness or 0s in the source marginals. Using uniform distribution
WARNING: Coupling at `(1.5, 2.0)` contains `1` empty row(s), e.g., due to unbalancedness or 0s in the source marginals. Using uniform distribution
ValueError: Transition matrix is not row stochastic, 5 rows do not sum to 1.
```

## 6

The repair carried into the next version (`fibdec_pseudotime_v2.ipynb`, same objects), wrapped around the CellRank call; it fires on 4 rows.

```bash
git -C $NH/xmens_final show b9fa19c:notebooks/fibdec_pseudotime_v2.ipynb | cell_src - 30 8
git -C $NH/xmens_final show b9fa19c:notebooks/fibdec_pseudotime_v2.ipynb | cell_out - 30
```
```
# cellrank's transition_matrix setter rejects any matrix whose rows don't sum to 1
# (rtol=1e-12), and the stitched moscot couplings routinely leave a few all-zero
# rows — cells that no coupling reached. This is the workaround from
# realtime_informed_pseudotime_human.ipynb cells 28-30, which is LOAD-BEARING, not
# diagnostic: catch the matrix before the setter, repair dead rows with a
# within-timepoint uniform, renormalize exactly, then assign.
#
# Mouse (5 transitions) hit this with 6 dead rows. Human (1 transition) did not.
Using automatic `threshold=3.0828566215145976e-44`
RealTimeKernel rejected the matrix: Transition matrix is not row stochastic, 4 rows do not sum to 1.
applying the documented repair ...
  dead rows (zero / non-finite): 4
  after repair: row sums [0.999999999999999, 1.000000000000000]
  non-stochastic (rtol=1e-12): 0
  assigned successfully

transition matrix (43197, 43197)  row sums [1.000000000, 1.000000000]
```

## 7

The repair in the xmens_final pipeline script, which logs the number of repaired rows.

```bash
git -C $NH/xmens_final grep -n 'row stochastic\|repaired rows' HEAD -- workflow/scripts/pseudotime/harmonized_pseudotime.py
```
```
HEAD:workflow/scripts/pseudotime/harmonized_pseudotime.py:145:    if "row stochastic" not in str(e):
HEAD:workflow/scripts/pseudotime/harmonized_pseudotime.py:168:    f"earlier {mass['earlier']:.3f}; repaired rows {n_dead}")
```

## 8

The human run of record on the 2025 objects (30 Oct 2026) repaired no rows. (No repaired-rows count was found for the mouse run of record on the 2025 objects, 22 Oct 2026; step 9 has the Nov 6 re-run on that object.)

```bash
grep -o '"repaired_rows": *[0-9]*' $NH/xmens_final/data/interim/human_pt_2026_10_30_run4_a40/provenance.json
```
```
"repaired_rows": 0
```

## 9

The 6 Nov 2026 re-runs on the 2025 objects (human 47,533 twice, mouse 42,290), with the repair still in the code: 0 rows each.

```bash
git -C $NH/xmens grep -h -o 'transition matrix.*repaired rows [0-9]*' eada0a5 -- notebooks/pseudotime_runs/
```
```
transition matrix (47533, 47533); mass later 0.395 / same 0.605 / earlier 0.000; repaired rows 0
transition matrix (47533, 47533); mass later 0.395 / same 0.605 / earlier 0.000; repaired rows 0
transition matrix (42290, 42290); mass later 0.697 / same 0.303 / earlier 0.000; repaired rows 0
```

## 10

The one November run where it fired: an interactive notebook last saved 7 Nov 2026, 18:39, on a copy of the November mouse object (41,811 nuclei), writing to a netscratch test path, not to an object of record. 5 rows repaired.

```bash
ls -l --time-style=+%F_%H:%M $NH/xmens_scripts/notebooks/pseudotime_mouse.ipynb | awk '{print $6, $7}'
cell_src $NH/xmens_scripts/notebooks/pseudotime_mouse.ipynb 1 4
cell_out $NH/xmens_scripts/notebooks/pseudotime_mouse.ipynb 32
```
```
2026-11-07_18:39 /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks/pseudotime_mouse.ipynb
# Parameters (papermill can override these)
INPUT_OBJECT = "/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad"
INPUT_OBJECT = "mouse_fibs_final.h5ad"
OUTPUT_OBJECT = "/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_pseudotime_notebooks/Cevrim_XMens_Only_Fib_Dec_Cells_pseudotime.h5ad"
5 rows kept same-day transitions only; 0 set uniform within their timepoint
transition matrix (41811, 41811); mass later 0.696 / same 0.304 / earlier 0.000; repaired rows 5
```

## 11

The runs that wrote the objects of record have no repair: the comment that replaced it, then each run's transition-matrix line. CellRank accepted every matrix (a rejected one would have stopped the run).

```bash
cell_src $NH/xmens/pseudotime/pseudotime_mouse_final_1110.ipynb 14 40 | grep -A3 'CellRank requires'
cell_out $NH/xmens/pseudotime/pseudotime_mouse_final_1106.ipynb 15 | grep 'transition matrix'
cell_out $NH/xmens/pseudotime/pseudotime_human_final_1107.ipynb 15 | grep 'transition matrix'
cell_out $NH/xmens/pseudotime/pseudotime_mouse_final_1110.ipynb 14 | grep 'transition matrix'
```
```
# CellRank requires every row to sum to 1 and raises "not row stochastic" otherwise (an underflowed transport row
# would do it). No run of record has hit this (repaired rows 0 in all), so nothing is repaired here: a run that
# hits it stops, and the affected nuclei get looked at before anything is changed.
rtk = cr.kernels.RealTimeKernel.from_moscot(tp)
transition matrix (41811, 41811); mass later 0.696 / same 0.304 / earlier 0.000
transition matrix (44949, 44949); mass later 0.409 / same 0.591 / earlier 0.000
transition matrix (41791, 41791); mass later 0.696 / same 0.304 / earlier 0.000
```

## 12

Those are the runs that wrote the three `_nick_pdt` objects (time and nuclei stored in each object).

```bash
$NH/.uv/envs/rsc_uv/bin/python -c "
import h5py
for n in ('20261106_mouse_fibdec_final_nick_pdt', '20261107_Human_fibdec_final_nick_pdt',
          '20261110_mouse_only_eFib_dec_final_nick_pdt'):
    g = h5py.File('$D/' + n + '.h5ad', 'r')['uns/pseudotime']
    print(n, g['computed'][()].decode(), int(g['n_cells'][()]), 'nuclei')
"
```
```
20261106_mouse_fibdec_final_nick_pdt 2026-11-10T10:03:35 41811 nuclei
20261107_Human_fibdec_final_nick_pdt 2026-11-10T11:51:47 44949 nuclei
20261110_mouse_only_eFib_dec_final_nick_pdt 2026-11-10T11:52:14 41791 nuclei
```

## 13

The settings were the same in the failing runs (step 5) and the runs of record: what changed is the input (which nuclei, and the scVI latent fit on them), not the parameters.

```bash
$NH/.uv/envs/rsc_uv/bin/python -c "
import h5py
g = h5py.File('$D/20261110_mouse_only_eFib_dec_final_nick_pdt.h5ad', 'r')['uns/pseudotime/params']
print({k: g[k][()].decode() for k in ('epsilon', 'tau', 'scale_cost', 'conn_weight', 'self_transitions', 'threshold')})
"
```
```
{'epsilon': '0.001', 'tau': '0.95', 'scale_cost': 'mean', 'conn_weight': '0.2', 'self_transitions': 'all', 'threshold': 'auto'}
```
