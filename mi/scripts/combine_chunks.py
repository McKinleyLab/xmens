"""Step 3: one species' per-chunk results -> one table with Benjamini-Hochberg q.

Reads the --n-chunks files mi_<chunk>.csv in --chunks and stops if any is missing; keeps the first row of each gene;
adds q, the Benjamini-Hochberg adjusted p over all genes (statsmodels fdr_bh); sorts by q. Writes --out, a CSV with
columns gene, mi (nats), p, q.
Same code as `mi_pseudotime.py merge` of the run of record (xmens_final/workflow/scripts/mi/), except that it stops when
a chunk is missing (the record's merge warned and went on).
"""
import argparse
import glob
import os
import sys

import pandas as pd
from statsmodels.stats.multitest import multipletests

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--chunks", required=True, help="folder of step 2's mi_<chunk>.csv files for one species")
parser.add_argument("--n-chunks", type=int, required=True, help="number of chunks expected")
parser.add_argument("--fdr", type=float, default=0.05, help="q threshold for the count printed at the end (default: 0.05)")
parser.add_argument("--out", required=True, help="output CSV; must not exist yet")
args = parser.parse_args()

if os.path.exists(args.out):
    sys.exit(f"{args.out} exists; not overwriting it")

# Files in name order, as in the run of record (mi_100.csv sorts before mi_1000.csv); the order of genes with equal q
# in the output follows it.
files = sorted(glob.glob(os.path.join(args.chunks, "mi_*.csv")))
expected = {f"mi_{chunk:03d}.csv" for chunk in range(args.n_chunks)}
found = {os.path.basename(path) for path in files}
if found != expected:
    missing, extra = sorted(expected - found), sorted(found - expected)
    sys.exit(f"{args.chunks}: {len(missing)} of {args.n_chunks} chunks missing {missing[:5]}, "
             f"{len(extra)} unexpected files {extra[:5]}")

table = pd.concat([pd.read_csv(path, usecols=["gene", "mi", "p"]) for path in files], ignore_index=True)
table = table.drop_duplicates("gene", keep="first")
table["q"] = multipletests(table["p"].to_numpy(), method="fdr_bh")[1]
table = table.sort_values("q").reset_index(drop=True)

# Written under a temporary name and renamed when complete, so an interrupted run never leaves a file that looks
# finished.
partial = os.path.join(os.path.dirname(args.out), ".partial." + os.path.basename(args.out))
table.to_csv(partial, index=False)
os.replace(partial, args.out)
print(f"{len(files)} chunks, {len(table):,} genes, {(table['q'] < args.fdr).sum():,} with q < {args.fdr}; "
      f"wrote {args.out}", flush=True)
