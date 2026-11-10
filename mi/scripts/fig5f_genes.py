"""Step 6: the one-to-one orthologs drawn in the cross-species pseudotime heatmap, Fig. 5F (Methods).

1. Pool: one-to-one pairs with both genes tested (9,638 in the run of record); per species, the top 10% of these pairs
   by MI (rank over the pairs); pool = pairs in the top 10% in both species (273).
2. Per species and pool gene: log1p(CP10K) per nucleus from raw counts (total counts over all genes), mean in 50
   equal-count bins of the pseudotime (pd.qcut on rank(method="first")), z-scored per gene (population SD);
   r = Pearson correlation of the human and mouse z-profiles (also over 20 bins, reported as profile_r_20bins).
3. Centroid per species = sum_b max(z_b, 0) * b / 49 / sum_b max(z_b, 0); gene centroid = mean of the two species;
   phase early (< 0.4), late (> 0.6) or mid.
4. Keep r > 0.8 (92 genes); take the 14 early and 14 late genes with the highest r; order them by centroid.
Writes the folder --out-dir: conserved_profiles_selected_genes.csv (the 28 genes, in drawing order),
conserved_profiles_pool_concordance.csv (every pool gene), selection_summary.json, and percell_mouse_pool.tsv.gz /
percell_human_pool.tsv.gz (raw counts of the pool genes per nucleus with a pseudotime, the input of the Monocle smoothing
that draws the heatmap).
Same code as crossspecies_gene_selection.py of the run of record (xmens_final/workflow/scripts/mi/), with the MI tables
and objects given as files and obs['pseudotime'] of the final objects (the nuclei and pseudotime the MI scored); the
comparison with the 28 genes drawn in the figure is made by compare_with_record.
"""
import argparse
import json
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

N_BINS_SELECTION = 50   # bins of the profile correlation used for the selection
N_BINS_CHECK = 20       # coarser bins, reported as profile_r_20bins
CP10K = 1e4

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--mi-mouse", required=True, help="step 3's mi_mouse.csv")
parser.add_argument("--mi-human", required=True, help="step 3's mi_human.csv")
parser.add_argument("--pairs", required=True, help="step 4's ortholog_pairs_one2one.csv")
parser.add_argument("--h5ad-mouse", required=True, help="final mouse object (raw counts in layers['counts'])")
parser.add_argument("--h5ad-human", required=True, help="final human object")
parser.add_argument("--pt-key", default="pseudotime", help="obs column holding the pseudotime (default: pseudotime)")
parser.add_argument("--top-fraction", type=float, default=0.10, help="top fraction by MI in both species (default: 0.10)")
parser.add_argument("--rmin", type=float, default=0.8, help="keep genes with profile r above this (default: 0.8)")
parser.add_argument("--n-early", type=int, default=14, help="early genes drawn (default: 14)")
parser.add_argument("--n-late", type=int, default=14, help="late genes drawn (default: 14)")
parser.add_argument("--out-dir", required=True, help="output folder; must not exist yet")
args = parser.parse_args()

out_dir = args.out_dir.rstrip("/")
if os.path.exists(out_dir):
    sys.exit(f"{out_dir} exists; not overwriting it")


def conserved_pool(mi_tables, pairs, top_fraction):
    """Pairs measured in both species, and the pool of pairs in the top fraction of MI in both."""
    mi = {species: pd.read_csv(path).drop_duplicates("gene").set_index("gene")["mi"]
          for species, path in mi_tables.items()}
    pairs = pairs[pairs["human"].isin(mi["human"].index) & pairs["mouse"].isin(mi["mouse"].index)]
    pairs = pairs.reset_index(drop=True)
    mi_mouse = mi["mouse"].reindex(pairs["mouse"]).to_numpy()
    mi_human = mi["human"].reindex(pairs["human"]).to_numpy()
    is_finite = np.isfinite(mi_mouse) & np.isfinite(mi_human)
    pairs = pairs[is_finite].reset_index(drop=True)
    mi_mouse, mi_human = mi_mouse[is_finite], mi_human[is_finite]

    k_top = int(round(top_fraction * len(pairs)))
    in_top = (
        (pd.Series(mi_mouse).rank(ascending=False).to_numpy() <= k_top)
        & (pd.Series(mi_human).rank(ascending=False).to_numpy() <= k_top)
    )
    pool = pairs[in_top].assign(mi_mouse=mi_mouse[in_top], mi_human=mi_human[in_top]).reset_index(drop=True)
    return pairs, k_top, pool


def per_nucleus_counts(h5ad, pt_key, genes):
    """Total counts, pseudotime and the raw counts of the pool genes, one row per nucleus with a pseudotime.

    Nuclei without a pseudotime (mouse apoptotic decidual nuclei) are dropped, as in the MI; rows keep object order.
    """
    adata = ad.read_h5ad(h5ad, backed="r")
    pt = adata.obs[pt_key].to_numpy(dtype=float)
    has_pt = np.isfinite(pt)
    columns = pd.Index(adata.var_names.astype(str)).get_indexer(pd.Index(genes))
    if (columns < 0).any():
        raise ValueError(f"{h5ad}: {int((columns < 0).sum())} pool genes missing from var")

    counts = sp.csr_matrix(adata.layers["counts"])[np.flatnonzero(has_pt)]
    total = np.asarray(counts.sum(axis=1)).ravel()
    table = pd.DataFrame(
        counts[:, columns].toarray(),
        columns=[f"cnt_{gene}" for gene in genes],
        index=pd.Index(adata.obs_names[has_pt].astype(str), name="cell"),
    )
    adata.file.close()
    table.insert(0, "total_counts", total)
    table.insert(1, "pt_rank_norm", pt[has_pt])
    print(f"{h5ad}: {len(table):,} nuclei with {pt_key} ({int((~has_pt).sum()):,} without)", flush=True)
    return table


def profiles(per_nucleus, genes, nbin):
    """Mean log1p(CP10K) per equal-count pseudotime bin (genes x bins) and its per-gene z-score."""
    counts = per_nucleus[[f"cnt_{gene}" for gene in genes]].to_numpy(float)
    expr = np.log1p(counts / per_nucleus["total_counts"].to_numpy(float)[:, None] * CP10K)
    bins = pd.qcut(per_nucleus["pt_rank_norm"].rank(method="first"), nbin, labels=False).to_numpy()
    means = np.vstack([expr[bins == b].mean(axis=0) for b in range(nbin)]).T
    sd = means.std(axis=1, ddof=0)
    sd[sd == 0] = 1
    return means, (means - means.mean(axis=1, keepdims=True)) / sd[:, None]


def centroid(z):
    """Pseudotime centroid on a 0-1 scale, weighted by the positive part of the z-scored profile."""
    weights = np.clip(z, 0, None)
    position = np.arange(z.shape[1]) / (z.shape[1] - 1)
    return (weights * position).sum(1) / weights.sum(1)


def rowwise_pearson(a, b):
    return np.array([np.corrcoef(a[i], b[i])[0, 1] for i in range(len(a))])


def concordance_table(pool, per_nucleus_mouse, per_nucleus_human):
    """Profile concordance, centroid, peak, range and phase of every pool pair."""
    human_genes, mouse_genes = list(pool["human"]), list(pool["mouse"])
    Mh50, Zh50 = profiles(per_nucleus_human, human_genes, N_BINS_SELECTION)
    Mm50, Zm50 = profiles(per_nucleus_mouse, mouse_genes, N_BINS_SELECTION)
    _, Zh20 = profiles(per_nucleus_human, human_genes, N_BINS_CHECK)
    _, Zm20 = profiles(per_nucleus_mouse, mouse_genes, N_BINS_CHECK)

    table = pool.copy()
    table["r50"] = rowwise_pearson(Zh50, Zm50)
    table["r20"] = rowwise_pearson(Zh20, Zm20)
    table["centroid_human"] = centroid(Zh50)
    table["centroid_mouse"] = centroid(Zm50)
    table["pt_centroid"] = (table["centroid_human"] + table["centroid_mouse"]) / 2
    table["peak_h"] = Zh50.argmax(1)
    table["peak_m"] = Zm50.argmax(1)
    table["range_h"] = Mh50.max(1) - Mh50.min(1)
    table["range_m"] = Mm50.max(1) - Mm50.min(1)
    # Phase cutoffs as written in the Methods: < 0.4 early, > 0.6 late.
    table["phase"] = np.where(table.pt_centroid < 0.4, "early", np.where(table.pt_centroid > 0.6, "late", "mid"))
    return table


pairs = pd.read_csv(args.pairs)[["human", "mouse"]]
universe, k_top, pool = conserved_pool({"mouse": args.mi_mouse, "human": args.mi_human}, pairs, args.top_fraction)
print(f"universe {len(universe):,} one-to-one pairs; top {args.top_fraction:.0%} = {k_top:,} per species; "
      f"pool {len(pool):,}", flush=True)

per_nucleus_mouse = per_nucleus_counts(args.h5ad_mouse, args.pt_key, list(pool["mouse"]))
per_nucleus_human = per_nucleus_counts(args.h5ad_human, args.pt_key, list(pool["human"]))
table = concordance_table(pool, per_nucleus_mouse, per_nucleus_human)

# Selection: the most concordant early and late genes above rmin, drawn in centroid order.
concordant = table[table["r50"] > args.rmin]
picked = pd.concat([
    concordant[concordant["phase"] == "early"].nlargest(args.n_early, "r50"),
    concordant[concordant["phase"] == "late"].nlargest(args.n_late, "r50"),
]).sort_values("pt_centroid")
selected = picked.rename(columns={"r50": "profile_r_50bins", "r20": "profile_r_20bins"})[
    ["human", "mouse", "profile_r_50bins", "profile_r_20bins", "phase", "pt_centroid", "peak_h", "peak_m"]
]
summary = {
    "mi_tables": {"mouse": args.mi_mouse, "human": args.mi_human},
    "universe": int(len(universe)),
    "k_top": k_top,
    "pool": int(len(pool)),
    "n_r_gt": int(len(concordant)),
    "rmin": args.rmin,
    "phase_counts_r_gt": concordant["phase"].value_counts().to_dict(),
    "median_r50_pool": round(float(table["r50"].median()), 3),
    "n_selected": int(len(selected)),
    "r50_range_selected": [round(float(selected["profile_r_50bins"].min()), 3),
                           round(float(selected["profile_r_50bins"].max()), 3)],
    "centroid_definition": "per species: sum_b max(z_b,0)*b/(nbin-1) / sum_b max(z_b,0) on 50-bin z-profiles; "
                           "gene = mean of species",
    "selected_human": list(selected["human"]),
}

# Written into a temporary folder that is renamed when complete, so an interrupted run never leaves a folder that looks
# finished (a leftover temporary folder is overwritten by the next run).
partial = os.path.join(os.path.dirname(out_dir), ".partial." + os.path.basename(out_dir))
os.makedirs(partial, exist_ok=True)
per_nucleus_mouse.to_csv(os.path.join(partial, "percell_mouse_pool.tsv.gz"), sep="\t")
per_nucleus_human.to_csv(os.path.join(partial, "percell_human_pool.tsv.gz"), sep="\t")
table.sort_values("r50", ascending=False).to_csv(os.path.join(partial, "conserved_profiles_pool_concordance.csv"),
                                                 index=False)
selected.to_csv(os.path.join(partial, "conserved_profiles_selected_genes.csv"), index=False)
with open(os.path.join(partial, "selection_summary.json"), "w") as f:
    json.dump(summary, f, indent=1)
os.rename(partial, out_dir)
print(json.dumps({key: value for key, value in summary.items() if key != "selected_human"}, indent=1))
print(f"selected ({len(selected)}, in drawing order): {', '.join(selected['human'])}", flush=True)
