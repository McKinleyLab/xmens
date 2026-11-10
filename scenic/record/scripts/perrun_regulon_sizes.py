#!/usr/bin/env python
"""Per-run vs consensus regulon size, both species.

For each of the 100 per-seed cisTarget outputs: how many regulons, and how many targets each has.
Compared against the consensus union at the floors of interest. Read-only; writes perrun_sizes.json/.txt.
"""
import glob, json, os, pickle
import numpy as np

LANE = "/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/scenic_v2_2026_08_28/scenic_seeds"
out = {}
for species in ("mouse", "human"):
    pkls = sorted(glob.glob(f"{LANE}/{species}/regulons_sketch15k/*.pkl"))
    assert len(pkls) == 100, f"{species}: {len(pkls)} pickles"
    n_reg, sizes, per_tf = [], [], {}
    for p in pkls:
        with open(p, "rb") as fh:
            regs = pickle.load(fh)
        seen = set()
        k = 0
        for r in regs:
            if r.name in seen:
                continue
            seen.add(r.name); k += 1
            s = len(r.gene2weight)
            sizes.append(s)
            per_tf.setdefault(r.name, []).append(s)
        n_reg.append(k)
    sizes = np.array(sizes)
    med_per_tf = {k: float(np.median(v)) for k, v in per_tf.items()}
    out[species] = {
        "n_pickles": len(pkls),
        "regulons_per_run": {"min": int(min(n_reg)), "median": float(np.median(n_reg)), "max": int(max(n_reg))},
        "target_count_per_regulon_all_runs": {q: float(np.quantile(sizes, q/100)) for q in (10, 25, 50, 75, 90)},
        "target_count_mean": float(sizes.mean()), "target_count_max": int(sizes.max()),
        "n_distinct_tfs_any_run": len(per_tf),
        "median_of_per_tf_median_size": float(np.median(list(med_per_tf.values()))),
        "per_tf_median_size": med_per_tf,
    }
    print(species, json.dumps({k: v for k, v in out[species].items() if k != "per_tf_median_size"}, indent=1), flush=True)

json.dump(out, open("perrun_sizes.json", "w"), indent=1)
with open("perrun_sizes.txt", "w") as f:
    for species, r in out.items():
        f.write(f"\n== {species}: {r['n_pickles']} runs, {r['regulons_per_run']['median']:.0f} regulons per run "
                f"(range {r['regulons_per_run']['min']}-{r['regulons_per_run']['max']}), "
                f"{r['n_distinct_tfs_any_run']} distinct TFs across runs\n")
        q = r["target_count_per_regulon_all_runs"]
        f.write(f"   targets per regulon per run: median {q[50]:.0f} (q10 {q[10]:.0f}, q25 {q[25]:.0f}, "
                f"q75 {q[75]:.0f}, q90 {q[90]:.0f}, max {r['target_count_max']}, mean {r['target_count_mean']:.0f})\n")
        f.write(f"   median across TFs of that TF's median per-run size: {r['median_of_per_tf_median_size']:.0f}\n")
print(open("perrun_sizes.txt").read())
