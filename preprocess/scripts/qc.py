"""Step 3b: quality control of one sample, matching the project QC notebook (= Cevrim et al. 2025).

Input: the step-3a AnnData (X = CellBender-corrected counts; obs has species/timepoint/coords; all cells present).

Procedure (Scanpy + DoubletDetection), faithful to the QC notebook:
  1. Compute QC metrics (mito/ribo/hb flags; mito on both mt-/MT- prefixes).
  2. Hard floor: >= min_counts (500) total counts.
  3. Mito %: fixed upper bound (8%).
  4. MAD outliers on total_counts and n_genes_by_counts, thresholded on a
     LOG10(x+1) scale with raw MAD, nmads=3.
  5. DoubletDetection (BoostClassifier, leiden, standard_scaling, 40 iters,
     n_top_var_genes=3000, p_thresh=1e-7, voter_thresh=0.5); remove doublets.
     Runs in the DoubletDetection 4.2 environment and stops if the package is missing or the fit fails.

Writes one QC'd .h5ad. Corrected counts kept in layers['counts']. Same code as
xmens_final/workflow/scripts/ingest_qc/qc.py (commit 59e00dd), with the Snakemake inputs replaced by the
command-line arguments below; their defaults are the settings of the run of record (record/config/config.yaml).
"""
import argparse
import os
import sys

import numpy as np
import scanpy as sc
from scipy.stats import median_abs_deviation


def number(text):
    """An int when the value is whole, so the log reads as the run of record's (e.g. "3-MAD", not "3.0-MAD")."""
    value = float(text)
    return int(value) if value.is_integer() else value


parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--in", dest="inp", required=True, help="step-3a .h5ad")
parser.add_argument("--out", required=True, help="output .h5ad; must not exist yet")
parser.add_argument("--sample", required=True, help="sample_id in samples.tsv")
parser.add_argument("--species", required=True, choices=["mouse", "human"])
parser.add_argument("--log", required=True)
parser.add_argument("--min-counts", type=number, default=500, help="config qc.min_counts")
parser.add_argument("--n-mads", type=number, default=3, help="config qc.n_mads")
parser.add_argument("--max-pct-mito", type=number, default=8, help="config qc.max_pct_mito")
parser.add_argument("--seed", type=int, default=0, help="config seed")
parser.add_argument("--n-iters", type=int, default=40, help="config doublet.n_iters")
parser.add_argument("--n-top-var-genes", type=int, default=3000, help="config doublet.n_top_var_genes")
parser.add_argument("--p-thresh", type=float, default=1e-7, help="config doublet.p_thresh")
parser.add_argument("--voter-thresh", type=float, default=0.5, help="config doublet.voter_thresh")
parser.add_argument("--n-jobs", type=int, default=-1, help="DoubletDetection threads; -1 (all) as in the record")
args = parser.parse_args()

inp = args.inp
out = args.out
sample = args.sample
species = args.species
logfile = args.log

qc = {"min_counts": args.min_counts, "n_mads": args.n_mads, "max_pct_mito": args.max_pct_mito}
seed = args.seed
dd_cfg = {"n_iters": args.n_iters, "n_top_var_genes": args.n_top_var_genes,
          "p_thresh": args.p_thresh, "voter_thresh": args.voter_thresh}

if os.path.exists(out):
    sys.exit(f"{out} exists; not overwriting it")

# Metrics MAD-thresholded on a log10(x+1) scale (matches notebook).
LOG_METRICS = ("total_counts", "n_genes_by_counts")


def _log10p(x):
    return np.log10(x + 1)


def _inv_log10p(y):
    return np.power(10, y) - 1


def outlier_bounds(values, nmads):
    """Raw-MAD bounds at +/- nmads from the median."""
    med = np.nanmedian(values)
    mad = median_abs_deviation(values, nan_policy="omit")  # scale=1.0 (raw)
    return med - nmads * mad, med + nmads * mad, mad


with open(logfile, "w") as log:
    def say(msg):
        print(f"[qc:{sample}] {msg}", file=log, flush=True)

    np.random.seed(seed)

    adata = sc.read_h5ad(inp)
    n0 = adata.n_obs
    say(f"loaded {n0} cells x {adata.n_vars} genes (species={species})")

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # --- 1. QC metrics (mito on both prefixes, like the notebook) ---
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
    adata.var["ribo"] = adata.var_names.str.startswith((
        "rps", "rpl", "mrps", "mrpl",
        "Rps", "Rpl", "Mrps", "Mrpl",
        "RPS", "RPL", "MRPS", "MRPL",
    ))
    adata.var["hb"] = (
        adata.var_names.str.contains("^Hb[^(p)]")
        | adata.var_names.str.contains("^HB[^(P)]")
    )
    n_mt = int(adata.var["mt"].sum())
    say(f"flagged {n_mt} mito genes (prefixes mt-/MT-)")
    if n_mt == 0:
        say("WARNING: no mito genes matched - check gene naming for this species")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"],
        percent_top=[20], log1p=True, inplace=True,
    )

    # Filtering is SEQUENTIAL: apply the
    # count floor first (so junk barcodes don't distort the MAD bounds), then
    # evaluate mito / MAD outliers on the surviving cells, in the notebook's
    # order (mito -> total_counts -> n_genes_by_counts).

    # --- 2. hard floor on total counts (applied first, mutates the object) ---
    min_counts = qc["min_counts"]
    n_low = int((adata.obs["total_counts"] < min_counts).sum())
    sc.pp.filter_cells(adata, min_counts=min_counts)
    say(f"min {min_counts} counts: removed {n_low}, {adata.n_obs}/{n0} remain")

    # --- 3. MAD outlier flags, computed on the FLOORED object ---
    # (bounds reflect the real cell population, not the junk below the floor)
    n_mads = qc["n_mads"]
    for metric in LOG_METRICS:
        vals = _log10p(adata.obs[metric].to_numpy())
        lo, hi, mad = outlier_bounds(vals, n_mads)
        if mad == 0:
            say(f"  WARNING: MAD is zero for {metric}; bounds collapse to the median")
        adata.obs[f"{metric}_outlier"] = (vals < lo) | (vals > hi)
        say(f"{n_mads}-MAD {metric} (log10): "
            f"bounds=[{_inv_log10p(lo):.1f}, {_inv_log10p(hi):.1f}], "
            f"flag {int(adata.obs[f'{metric}_outlier'].sum())}")

    # mito outlier flag (fixed upper bound)
    max_mito = qc["max_pct_mito"]
    adata.obs["pct_counts_mt_outlier"] = adata.obs["pct_counts_mt"] > max_mito
    say(f"mito > {max_mito}%: flag {int(adata.obs['pct_counts_mt_outlier'].sum())}")

    # --- 4. apply the flags sequentially, in notebook order ---
    for metric in ("pct_counts_mt", "total_counts", "n_genes_by_counts"):
        n_before = adata.n_obs
        adata = adata[~adata.obs[f"{metric}_outlier"]].copy()
        say(f"  filter {metric}: removed {n_before - adata.n_obs}, {adata.n_obs} remain")

    say(f"after threshold QC: {adata.n_obs}/{n0} cells")

    # --- 5. DoubletDetection ---
    # Imported after the seed, in the same order as the run of record; a missing package stops the run here.
    import doubletdetection
    ver = getattr(doubletdetection, "__version__", "?")
    n_iters = dd_cfg["n_iters"]
    n_top = dd_cfg["n_top_var_genes"]
    p_thresh = dd_cfg["p_thresh"]
    voter_thresh = dd_cfg["voter_thresh"]
    say(f"DoubletDetection {ver}: {n_iters} iters, n_top_var_genes={n_top}, "
        f"voter_thresh={voter_thresh}")
    clf = doubletdetection.BoostClassifier(
        n_iters=n_iters,
        clustering_algorithm="leiden",
        standard_scaling=True,
        n_top_var_genes=n_top,
        n_jobs=args.n_jobs,
    )
    # Fit on counts (matches notebook: clf.fit(adata.X) on corrected counts).
    preds = clf.fit(adata.layers["counts"]).predict(
        p_thresh=p_thresh, voter_thresh=voter_thresh
    )
    adata.obs["doubletdetection_predicted_doublet"] = preds.astype(bool)
    adata.obs["doubletdetection_doublet_score"] = clf.doublet_score()
    n_doub = int(adata.obs["doubletdetection_predicted_doublet"].sum())
    say(f"flagged {n_doub} doublets; removing")
    adata = adata[~adata.obs["doubletdetection_predicted_doublet"].to_numpy()].copy()

    say(f"final: {adata.n_obs} cells x {adata.n_vars} genes "
        f"({n0 - adata.n_obs} removed total)")
    # Written under a temporary name and renamed when complete, so an interrupted run never leaves a file that
    # looks finished.
    partial = os.path.join(os.path.dirname(out), "." + os.path.basename(out) + ".partial.h5ad")
    adata.write_h5ad(partial)
    os.replace(partial, out)
    say(f"wrote -> {out}")
