"""Step 2: mutual information (MI) between expression and pseudotime, with a permutation p-value, for one chunk of genes.

The genes of step 1 are split into --n-chunks consecutive chunks of near-equal size; this scores chunk --chunk.
  MI    scikit-learn's mutual_info_regression (k-nearest-neighbour estimator, 3 neighbours, random_state 0) on all
        nuclei, in nats
  null  --permutations shuffles of the pseudotime, drawn once for the chunk by numpy's default_rng(chunk) and shared by
        all of its genes; the MI of each shuffle is computed the same way
  p     (number of shuffles with MI >= observed + 1) / (permutations + 1)
Writes --out, a CSV with columns gene, mi, p.
Same code as `mi_pseudotime.py worker` of the run of record (xmens_final/workflow/scripts/mi/, run with --perms 1001,
which counts the observed pseudotime, i.e. 1,000 shuffles; --subsample-cells 0; seed 0), with the inputs read from
step 1's folder.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.feature_selection import mutual_info_regression

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--inputs", required=True, help="step 1's folder for one species")
parser.add_argument("--n-chunks", type=int, required=True, help="number of chunks the genes are split into")
parser.add_argument("--chunk", type=int, required=True, help="the chunk to score, 0 to n-chunks - 1")
parser.add_argument("--permutations", type=int, default=1000, help="shuffles of the pseudotime (default: 1000)")
parser.add_argument("--seed", type=int, default=0,
                    help="MI estimator seed; the shuffles are seeded with seed * 100,000 + chunk (default: 0)")
parser.add_argument("--out", required=True, help="output CSV; must not exist yet")
args = parser.parse_args()

if os.path.exists(args.out):
    sys.exit(f"{args.out} exists; not overwriting it")
if not 0 <= args.chunk < args.n_chunks:
    sys.exit(f"--chunk {args.chunk} is outside 0 to {args.n_chunks - 1}")

expression = np.load(os.path.join(args.inputs, "expression.npz"))["X"]           # nuclei x genes, float32
pseudotime = np.load(os.path.join(args.inputs, "pseudotime.npy")).astype(float)
with open(os.path.join(args.inputs, "genes.txt")) as f:
    genes = f.read().splitlines()
if not expression.shape == (len(pseudotime), len(genes)):
    sys.exit(f"inputs disagree: expression {expression.shape}, {len(pseudotime)} pseudotimes, {len(genes)} genes")

# Chunk c holds genes edges[c] to edges[c + 1] - 1.
edges = np.linspace(0, len(genes), args.n_chunks + 1, dtype=int)
first, end = int(edges[args.chunk]), int(edges[args.chunk + 1])
print(f"chunk {args.chunk} of {args.n_chunks}: genes {first} to {end - 1} ({end - first} genes), "
      f"{len(pseudotime):,} nuclei, {args.permutations} permutations", flush=True)

# The shuffled pseudotimes, drawn once and shared by every gene of the chunk; as in the run of record, each one shuffles
# the previous one in place.
rng = np.random.default_rng(args.seed * 100_000 + args.chunk)
shuffles = np.empty((args.permutations, len(pseudotime)))
shuffled = pseudotime.copy()
for i in range(args.permutations):
    rng.shuffle(shuffled)
    shuffles[i] = shuffled

rows = []
for g in range(first, end):
    x = expression[:, [g]]
    mi = float(mutual_info_regression(x, pseudotime, discrete_features=False, random_state=args.seed)[0])
    # The null MIs are kept as float32 and compared with the Python float mi, as in the run of record; keep both types,
    # or a tie at float32 precision can change p by 1/1,001.
    null = np.empty(args.permutations, dtype=np.float32)
    for i in range(args.permutations):
        null[i] = mutual_info_regression(x, shuffles[i], discrete_features=False, random_state=args.seed)[0]
    p = (float((null >= mi).sum()) + 1.0) / (args.permutations + 1.0)
    rows.append({"gene": genes[g], "mi": mi, "p": p})
    print(f"{genes[g]}: MI {mi:.4f} nats, p {p:.4g}", flush=True)

# Written under a temporary name and renamed when complete, so an interrupted task never leaves a file that looks
# finished.
partial = os.path.join(os.path.dirname(args.out), ".partial." + os.path.basename(args.out))
pd.DataFrame(rows, columns=["gene", "mi", "p"]).to_csv(partial, index=False)
os.replace(partial, args.out)
print(f"wrote {args.out}", flush=True)
