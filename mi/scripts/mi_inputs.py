"""Step 1: one species' final object -> the expression matrix, pseudotime and gene list that step 2 scores.

Keeps the nuclei with a pseudotime (in mouse, the 907 apoptotic decidual nuclei have none) and the genes detected
(count > 0) in at least 1% of those nuclei and with at least 10 counts in total, then computes log1p(CP10K) from the raw
counts of the kept genes: each nucleus is scaled to 10,000 counts summed over the kept genes, then log(1 + x).
Writes the folder --out-dir:
  expression.npz   X, nuclei x genes, float32
  pseudotime.npy   the pseudotime of each nucleus, float64
  genes.txt        gene names, one per line, in column order
Same computation as `mi_pseudotime.py prep` of the run of record (xmens_final/workflow/scripts/mi/, run with
--drop-missing-pt --min-detect-frac 0.01 --min-counts 10 and log1p(CP10K) derived from the counts), with the Snakemake
inputs replaced by the command-line arguments below and the object's own obs['pseudotime'].
"""
import argparse
import os
import sys

import anndata as ad
import numpy as np
import scipy.sparse as sp

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--h5ad", required=True, help="final object, with raw counts in layers['counts']")
parser.add_argument("--pt-key", default="pseudotime", help="obs column holding the pseudotime (default: pseudotime)")
parser.add_argument("--min-detect-frac", type=float, default=0.01,
                    help="keep genes detected in at least this fraction of nuclei (default: 0.01)")
parser.add_argument("--min-counts", type=int, default=10,
                    help="and with at least this many counts in total (default: 10)")
parser.add_argument("--out-dir", required=True, help="output folder; must not exist yet")
args = parser.parse_args()

out_dir = args.out_dir.rstrip("/")
if os.path.exists(out_dir):
    sys.exit(f"{out_dir} exists; not overwriting it")

adata = ad.read_h5ad(args.h5ad)
print(f"{args.h5ad}: {adata.n_obs:,} nuclei x {adata.n_vars:,} genes", flush=True)

# Nuclei with a pseudotime.
pseudotime = adata.obs[args.pt_key].to_numpy(dtype=float)
has_pt = np.isfinite(pseudotime)
print(f"{has_pt.sum():,} nuclei with a pseudotime; {(~has_pt).sum():,} without, dropped", flush=True)
counts = sp.csr_matrix(adata.layers["counts"])[np.flatnonzero(has_pt)]
pseudotime = pseudotime[has_pt]
genes = adata.var_names.to_numpy().astype(str)

# Genes detected in at least 1% of these nuclei, then, among them, genes with at least 10 counts.
detected = np.asarray((counts > 0).mean(axis=0)).ravel() >= args.min_detect_frac
counts, genes = counts[:, detected], genes[detected]
print(f"{detected.sum():,} genes detected in at least {args.min_detect_frac:.0%} of nuclei", flush=True)
enough = np.asarray(counts.sum(axis=0)).ravel() >= args.min_counts
counts, genes = counts[:, enough], genes[enough]
print(f"{enough.sum():,} of them with at least {args.min_counts} counts: the genes tested", flush=True)

# log1p(CP10K) over the kept genes, in float64, then stored as float32 (the operations and their order are those of the
# run of record, so the matrix is identical to its expr.npz).
counts = counts.astype(np.float64)
library_size = np.asarray(counts.sum(axis=1)).ravel()
library_size[library_size == 0] = 1.0
expression = (sp.diags(1e4 / library_size) @ counts).toarray()
np.log1p(expression, out=expression)
expression = np.ascontiguousarray(expression, dtype=np.float32)

# Written into a temporary folder that is renamed when complete, so an interrupted run never leaves a folder that looks
# finished (a leftover temporary folder is overwritten by the next run).
partial = os.path.join(os.path.dirname(out_dir), ".partial." + os.path.basename(out_dir))
os.makedirs(partial, exist_ok=True)
np.savez_compressed(os.path.join(partial, "expression.npz"), X=expression)
np.save(os.path.join(partial, "pseudotime.npy"), pseudotime)
with open(os.path.join(partial, "genes.txt"), "w") as f:
    f.write("\n".join(genes) + "\n")
os.rename(partial, out_dir)
print(f"wrote {out_dir}: {expression.shape[0]:,} nuclei x {expression.shape[1]:,} genes", flush=True)
