#!/usr/bin/env python3
"""
mi_worker.py — run mutual information for one gene chunk. For use with
integrated pseudotime outputs, human or mouse.

Inputs
------
--h5ad         Path to .h5ad with `layers["log1p"]` and obs[--pt-key]
--chunk        Zero-based chunk index (defaults to $SLURM_ARRAY_TASK_ID)
--n-chunks     Total number of chunks; defaults to SLURM array size if present
--layer        Expression layer to use (default: log1p)
--pt-key       Pseudotime key in .obs (default: dpt_pseudotime)
--perms        Total permutations incl. observed (default: 1001)
--outdir       Output directory (default: ./results)

Output
------
CSV: {outdir}/mi_{chunk:03d}.csv with columns:
  gene, mi, p, mi_null_1, ..., mi_null_{perms-1}
TXT: {outdir}/chunk_{chunk:03d}.txt listing the genes in this chunk.
"""

import os, sys, math, argparse, pathlib
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm
from sklearn.feature_selection import mutual_info_regression

def empirical_p(obs, null):
    return ((null >= obs).sum() + 1) / (len(null) + 1)

def get_chunk_indices(total_genes, n_chunks, chunk_id):
    split = np.linspace(0, total_genes, n_chunks + 1, dtype=int)
    return split[chunk_id], split[chunk_id + 1]

def detect_slurm_chunk():
    return int(os.environ.get("SLURM_ARRAY_TASK_ID", "-1"))

def detect_slurm_nchunks():
    # Works on most SLURM setups: max - min + 1 (e.g., 0-999 → 1000)
    smin = os.environ.get("SLURM_ARRAY_TASK_MIN")
    smax = os.environ.get("SLURM_ARRAY_TASK_MAX")
    if smin is not None and smax is not None:
        return int(smax) - int(smin) + 1
    # Fallback: let user pass --n-chunks
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, type=pathlib.Path)
    ap.add_argument("--chunk", type=int, help="0-based chunk index")
    ap.add_argument("--n-chunks", type=int, help="total number of chunks")
    ap.add_argument("--layer", default="log1p")
    ap.add_argument("--pt-key", default="dpt_pseudotime")
    ap.add_argument("--perms", type=int, default=1001)
    ap.add_argument("--outdir", type=pathlib.Path, default="results")
    args = ap.parse_args()

    chunk = args.chunk if args.chunk is not None else detect_slurm_chunk()

    n_chunks = args.n_chunks if args.n_chunks is not None else detect_slurm_nchunks()

    print(f"[worker {chunk:03d}] reading {args.h5ad}")
    ad = sc.read_h5ad(args.h5ad)

    X = ad.layers[args.layer]


    # dense for sklearn MI
    X = X.toarray() if hasattr(X, "toarray") else X
    X = X.astype(np.float32, copy=False)

    y = ad.obs[args.pt_key].to_numpy()

    start, end = get_chunk_indices(ad.n_vars, n_chunks, chunk)
    genes = ad.var_names[start:end].to_numpy()
    Xc = X[:, start:end]

    print(f"[worker {chunk:03d}] {len(genes)} genes ({start}:{end})")

    rng = np.random.default_rng(seed=chunk)
    Ptot = int(args.perms)
    null_cols = [f"mi_null_{i}" for i in range(1, Ptot)]

    records = []
    for j, g in tqdm(list(enumerate(genes)), desc=f"chunk {chunk:03d}"):
        x = Xc[:, [j]]
        mi_obs = mutual_info_regression(x, y, discrete_features=False)[0]

        mi_null = np.empty(Ptot - 1, dtype=np.float32)
        y_shuf = y.copy()
        for p in range(Ptot - 1):
            rng.shuffle(y_shuf)
            mi_null[p] = mutual_info_regression(x, y_shuf, discrete_features=False)[0]

        records.append({"gene": g, "mi": float(mi_obs), "p": empirical_p(mi_obs, mi_null), **dict(zip(null_cols, mi_null))})

    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame.from_records(records, columns=["gene", "mi", "p"] + null_cols)\
      .to_csv(outdir / f"mi_{chunk:03d}.csv", index=False)
    np.savetxt(outdir / f"chunk_{chunk:03d}.txt", genes, fmt="%s")

    print(f"[worker {chunk:03d}] wrote {outdir}/mi_{chunk:03d}.csv")

if __name__ == "__main__":
    main()