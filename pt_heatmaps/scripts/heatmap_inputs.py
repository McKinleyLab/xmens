"""Step 1: the final objects, the MI tables and the gene lists -> the inputs of Monocle (step 2) and of the figures.

Per species, from the final object: the nuclei with a pseudotime (in mouse, the 907 apoptotic decidual nuclei have
none) and the raw counts of all their genes, which Monocle needs for its size factors and dispersion trend.
The genes drawn:
  Fig. 4H  the Yang et al. (2023) panel (--yang), in published order, matched to the mouse gene symbols
  Fig. 5F  the 28 genes of --fig5f-genes, in drawing order (xmens_mi step 6)
  MI       per species, the top pseudotime-associated genes among the one-to-one orthologs significant in both species:
           MI above the 90th percentile of the species' significant genes (q < 0.05); 1,133 mouse and 1,201 human
           genes, 373 pairs in both (the K, n and k of the Methods)
Writes the folder --out-dir:
  <species>/counts.mtx        raw counts, genes x nuclei (Matrix Market, integer)
  <species>/genes.txt         the genes, in the row order of counts.mtx
  cells_<species>.tsv.gz      the nuclei, in the column order of counts.mtx: cell, pt_rank_norm (= obs['pseudotime']),
                              label (cell_subtype), timepoint
  yang_panel.tsv              gene (published name), symbol (mouse gene), published_row
  fig5f_genes.csv             a copy of --fig5f-genes
  mi_top_<species>.csv        gene, partner (one-to-one ortholog), mi (nats), shared (top in both species); by MI
Same cells, counts and gene sets as export_heatmap_inputs_final_objects.py of the 2 Nov 2026 heatmap package
(record/), which wrote per-gene count tables and a dispersion table instead of the full counts.
"""
import argparse
import os
import shutil
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--h5ad-mouse", required=True)
parser.add_argument("--h5ad-human", required=True)
parser.add_argument("--mi-mouse", required=True, help="MI table of the mouse genes (gene, mi, p, q)")
parser.add_argument("--mi-human", required=True, help="MI table of the human genes")
parser.add_argument("--pairs", required=True, help="one-to-one ortholog pairs (columns human, mouse)")
parser.add_argument("--fig5f-genes", required=True, help="the Fig. 5F genes (columns human, mouse), in drawing order")
parser.add_argument("--yang", required=True, help="the Yang et al. (2023) panel (column gene), in published order")
parser.add_argument("--out-dir", required=True, help="output folder; must not exist yet")
args = parser.parse_args()

out_dir = args.out_dir.rstrip("/")
if os.path.exists(out_dir):
    sys.exit(f"{out_dir} exists; not overwriting it")
# Everything is written into a temporary folder that is renamed when complete, so an interrupted run never leaves a
# folder that looks finished (a leftover temporary folder is overwritten by the next run).
partial = os.path.join(os.path.dirname(out_dir), ".partial." + os.path.basename(out_dir))
os.makedirs(partial, exist_ok=True)

# ---- per species: the nuclei with a pseudotime and the raw counts of all genes ---------------------------------------
symbols = {}
for species, h5ad in (("mouse", args.h5ad_mouse), ("human", args.h5ad_human)):
    adata = ad.read_h5ad(h5ad, backed="r")
    pseudotime = adata.obs["pseudotime"].to_numpy(dtype=float)
    keep = np.isfinite(pseudotime)
    genes = adata.var_names.astype(str)
    if genes.duplicated().any():
        sys.exit(f"{h5ad}: gene names are not unique")
    counts = sp.csr_matrix(adata.layers["counts"])[np.flatnonzero(keep)]
    counts.data = np.round(counts.data).astype(np.int64)
    counts.sort_indices()
    obs = adata.obs.loc[keep]
    adata.file.close()
    print(f"{species}: {keep.sum():,} nuclei with a pseudotime ({(~keep).sum():,} without), {len(genes):,} genes, "
          f"{counts.nnz:,} non-zero counts", flush=True)

    os.makedirs(os.path.join(partial, species))
    scipy.io.mmwrite(os.path.join(partial, species, "counts.mtx"), counts.T.tocoo(), field="integer",
                     symmetry="general")
    with open(os.path.join(partial, species, "genes.txt"), "w") as f:
        f.write("\n".join(genes) + "\n")
    cells = pd.DataFrame({"cell": obs.index.astype(str), "pt_rank_norm": pseudotime[keep],
                          "label": obs["cell_subtype"].astype(str).to_numpy(),
                          "timepoint": obs["timepoint"].astype(str).to_numpy()})
    cells.to_csv(os.path.join(partial, f"cells_{species}.tsv.gz"), sep="\t", index=False, compression="gzip")
    symbols[species] = genes

# ---- Fig. 4H: the Yang panel, matched to the mouse gene symbols ------------------------------------------------------
yang = pd.read_csv(args.yang, sep="\t")
mouse_symbol = {gene.upper(): gene for gene in symbols["mouse"]}
missing = [gene for gene in yang["gene"] if gene.upper() not in mouse_symbol]
if missing:
    sys.exit(f"Yang panel genes missing from the mouse object: {missing}")
yang.insert(1, "symbol", [mouse_symbol[gene.upper()] for gene in yang["gene"]])
yang.to_csv(os.path.join(partial, "yang_panel.tsv"), sep="\t", index=False)

# ---- Fig. 5F: the 28 genes, in drawing order -------------------------------------------------------------------------
fig5f = pd.read_csv(args.fig5f_genes)
for species in ("mouse", "human"):
    missing = sorted(set(fig5f[species]) - set(symbols[species]))
    if missing:
        sys.exit(f"Fig. 5F genes missing from the {species} object: {missing}")
shutil.copyfile(args.fig5f_genes, os.path.join(partial, "fig5f_genes.csv"))

# ---- MI heatmaps: the top genes of each species among the one-to-one orthologs significant in both -------------------
mi = {"mouse": pd.read_csv(args.mi_mouse).drop_duplicates("gene").set_index("gene"),
      "human": pd.read_csv(args.mi_human).drop_duplicates("gene").set_index("gene")}
significant = {species: table[table["q"] < 0.05] for species, table in mi.items()}
pairs = pd.read_csv(args.pairs)
background = pairs[pairs["mouse"].isin(significant["mouse"].index)
                   & pairs["human"].isin(significant["human"].index)].copy()
background["mi_mouse"] = significant["mouse"]["mi"].reindex(background["mouse"]).to_numpy()
background["mi_human"] = significant["human"]["mi"].reindex(background["human"]).to_numpy()
threshold = {species: float(np.percentile(table["mi"], 90)) for species, table in significant.items()}
top_mouse = background["mi_mouse"] > threshold["mouse"]
top_human = background["mi_human"] > threshold["human"]
both = top_mouse & top_human
print(f"MI heatmaps: N {len(background):,} pairs, {top_mouse.sum():,} top mouse genes (K), {top_human.sum():,} top "
      f"human genes (n), {both.sum():,} in both (k)", flush=True)
for species, partner, top in (("mouse", "human", top_mouse), ("human", "mouse", top_human)):
    table = background[top].assign(gene=background[species], partner=background[partner],
                                    mi=background[f"mi_{species}"], shared=both[top])
    table = table[["gene", "partner", "mi", "shared"]].sort_values("mi", ascending=False)
    missing = sorted(set(table["gene"]) - set(symbols[species]))
    if missing:
        sys.exit(f"{len(missing)} MI-heatmap genes missing from the {species} object, e.g. {missing[:5]}")
    table.to_csv(os.path.join(partial, f"mi_top_{species}.csv"), index=False)

os.rename(partial, out_dir)
print(f"wrote {out_dir}", flush=True)
