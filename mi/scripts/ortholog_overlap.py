"""Step 4: overlap of the human and mouse top MI genes among one-to-one orthologs (Methods; Supplementary Fig. S29C-D).

1. One-to-one pairs: Ensembl Compara components with homology type ortholog_one2one and biotype protein_coding, one
   human and one mouse gene each; genes found in more than one pair are dropped.
2. Per species, the significant genes (q < 0.05) and the top genes: significant genes whose MI is above the 90th
   percentile of the significant genes' MI (numpy percentile, linear interpolation).
3. N: pairs with both genes significant. K and n: the pairs among them whose mouse or human gene is a top gene; k: the
   pairs top in both species. Expected overlap K n / N; p = P(X >= k) for X hypergeometric with parameters N, K, n,
   computed on the log scale (scipy hypergeom.logsf) so that it cannot underflow.
4. Steps 2 and 3 at every percentile cutoff from the 40th to the 99th, with the 2.5-97.5% range of the hypergeometric
   null, expressed as fold enrichment.
Writes to --out-dir: ortholog_pairs_one2one.csv, methods_overlap.json and overlap_percentile_sweep.csv.
Same code as the functions one_to_one_pairs, methods_overlap and percentile_sweep of mi_overlap.py of the run of record
(xmens_final/workflow/scripts/mi/). Its other outputs (control arms, shuffled-pairing null, 1-30% cutoff sweep,
shared-universe specificity) are not cited in the Methods and are not ported.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--mi-mouse", required=True, help="step 3's mi_mouse.csv")
parser.add_argument("--mi-human", required=True, help="step 3's mi_human.csv")
parser.add_argument("--orthologs", required=True, help="Ensembl Compara human-mouse components (TSV)")
parser.add_argument("--q-thresh", type=float, default=0.05, help="MI significance threshold on q (default: 0.05)")
parser.add_argument("--top-percentile", type=float, default=90.0,
                    help="top genes: MI above this percentile of the significant genes' MI (default: 90)")
parser.add_argument("--out-dir", required=True, help="output folder; must not exist yet")
args = parser.parse_args()

out_dir = args.out_dir.rstrip("/")
if os.path.exists(out_dir):
    sys.exit(f"{out_dir} exists; not overwriting it")


def one_to_one_pairs(components):
    """Strict one-to-one protein-coding (human, mouse) pairs, one row per Compara component."""
    needed = {"component_id", "species", "gene_symbol", "homology_type", "biotype"}
    if not needed <= set(components.columns):
        raise ValueError(f"ortholog table lacks {sorted(needed - set(components.columns))}")
    one2one = components[
        (components["homology_type"] == "ortholog_one2one") & (components["biotype"] == "protein_coding")
    ]
    if not (one2one.groupby(["component_id", "species"]).size() == 1).all():
        raise ValueError("a one-to-one component has more than one gene in a species")
    wide = one2one.pivot_table(index="component_id", columns="species", values="gene_symbol", aggfunc="first")
    pairs = wide.dropna()[["human", "mouse"]].drop_duplicates()
    is_unique = ~pairs["human"].duplicated(keep=False) & ~pairs["mouse"].duplicated(keep=False)
    return pairs[is_unique]


def overlap_at(pairs, mi_mouse, mi_human, percentile):
    """N, K, n, k, the expected overlap, fold and log10 p of the top sets at one percentile cutoff."""
    sets = {}
    for species, table in (("mouse", mi_mouse), ("human", mi_human)):
        significant = table[table["q"] < args.q_thresh]
        threshold = float(np.percentile(significant["mi"], percentile))
        sets[species] = {
            "tested": len(table),
            "significant": len(significant),
            "threshold_nats": threshold,
            "top": set(significant.index[significant["mi"] > threshold]),
            "significant_genes": set(significant.index),
        }
    background = pairs[
        pairs["mouse"].isin(sets["mouse"]["significant_genes"]) & pairs["human"].isin(sets["human"]["significant_genes"])
    ]
    is_top_mouse = background["mouse"].isin(sets["mouse"]["top"])
    is_top_human = background["human"].isin(sets["human"]["top"])
    N, K, n = len(background), int(is_top_mouse.sum()), int(is_top_human.sum())
    k = int((is_top_mouse & is_top_human).sum())
    expected = K * n / N
    log10_p = float(hypergeom.logsf(k - 1, N, K, n)) / np.log(10)
    return sets, N, K, n, k, expected, log10_p


mi_mouse = pd.read_csv(args.mi_mouse).drop_duplicates("gene").set_index("gene")
mi_human = pd.read_csv(args.mi_human).drop_duplicates("gene").set_index("gene")
pairs = one_to_one_pairs(pd.read_csv(args.orthologs, sep="\t"))
print(f"{len(pairs):,} one-to-one pairs from {args.orthologs}", flush=True)

# The Methods numbers, at the 90th percentile.
sets, N, K, n, k, expected, log10_p = overlap_at(pairs, mi_mouse, mi_human, args.top_percentile)
methods = {
    "q_thresh": args.q_thresh,
    "top_percentile": args.top_percentile,
    "N": N,
    "K": K,
    "n": n,
    "k": k,
    "expected": expected,
    "fold": k / expected,
    "log10_p": log10_p,
}
for species, values in sets.items():
    methods[species] = {
        "tested": values["tested"],
        "significant": values["significant"],
        "threshold_nats": values["threshold_nats"],
        "threshold_bits": values["threshold_nats"] / np.log(2),
        "top": len(values["top"]),
    }
methods["mi_tables"] = {"mouse": args.mi_mouse, "human": args.mi_human}
exponent = int(np.floor(log10_p))
print(f"N {N:,}, K {K:,}, n {n:,}, k {k:,}; expected {expected:.1f}, fold {k / expected:.2f}, "
      f"p = {10 ** (log10_p - exponent):.2f} x 10^{exponent} (log10 p {log10_p:.4f})", flush=True)

# The same at every percentile from the 40th to the 99th (Supplementary Fig. S29C-D).
sweep = []
for percentile in range(40, 100):
    sets, N, K, n, k, expected, log10_p = overlap_at(pairs, mi_mouse, mi_human, percentile)
    sweep.append({
        "percentile": percentile,
        "thr_mouse": sets["mouse"]["threshold_nats"],
        "thr_human": sets["human"]["threshold_nats"],
        "top_mouse_genomewide": len(sets["mouse"]["top"]),
        "top_human_genomewide": len(sets["human"]["top"]),
        "N": N,
        "K_mouse": K,
        "n_human": n,
        "overlap": k,
        "expected": expected,
        "fold": k / expected,
        "neglog10_p": -log10_p,
        "null_fold_lo": hypergeom.ppf(0.025, N, K, n) / expected,
        "null_fold_hi": hypergeom.ppf(0.975, N, K, n) / expected,
        "pct_of_human_top": 100 * k / n,
    })

# Written into a temporary folder that is renamed when complete, so an interrupted run never leaves a folder that looks
# finished (a leftover temporary folder is overwritten by the next run).
partial = os.path.join(os.path.dirname(out_dir), ".partial." + os.path.basename(out_dir))
os.makedirs(partial, exist_ok=True)
pairs.to_csv(os.path.join(partial, "ortholog_pairs_one2one.csv"), index=False)
with open(os.path.join(partial, "methods_overlap.json"), "w") as f:
    f.write(json.dumps(methods, indent=1) + "\n")
pd.DataFrame(sweep).to_csv(os.path.join(partial, "overlap_percentile_sweep.csv"), index=False)
os.rename(partial, out_dir)
print(f"wrote {out_dir}", flush=True)
