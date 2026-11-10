"""Step 6: QC figures for one sample, with the thresholds qc.py applied drawn on the distributions.

Reads the step-3a object (every CellBender cell, before QC) and the step-3b object (after QC) and writes one figure,
<sample>_qc.png, with four panels: total counts and detected genes per nucleus (log scale, with the 500-count floor and
the 3-MAD bounds), the mitochondrial percentage (with the 8% bound) and the DoubletDetection scores of the nuclei that
passed. The bounds are recomputed exactly as qc.py computes them: on the nuclei above the floor, from the raw MAD of
log10(x + 1). With --qc-log they are checked against the bounds qc.py wrote to its log, so the lines on the figure
are the thresholds that were applied. Also writes summary.txt (nuclei before and after each filter).

Ported from xmens_final/workflow/scripts/ingest_qc/qc_report.py (record/3_qc/qc_report.py), which was never run;
that version computed the bounds before the floor, which is not what qc.py does.
"""
import argparse
import os
import re
import sys

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from scipy.stats import median_abs_deviation

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--raw", required=True, help="step-3a .h5ad (every CellBender cell)")
parser.add_argument("--qc", required=True, help="step-3b .h5ad (nuclei that passed QC)")
parser.add_argument("--sample", required=True, help="sample_id in samples.tsv")
parser.add_argument("--species", required=True, choices=["mouse", "human"])
parser.add_argument("--outdir", required=True, help="folder for <sample>_qc.png and summary.txt; created if needed")
parser.add_argument("--log", required=True)
parser.add_argument("--qc-log", help="qc.py's log for this sample; the recomputed bounds must match its bounds lines")
parser.add_argument("--min-counts", type=float, default=500, help="as qc.py")
parser.add_argument("--n-mads", type=float, default=3, help="as qc.py")
parser.add_argument("--max-pct-mito", type=float, default=8, help="as qc.py")
args = parser.parse_args()

BLUE, RED, GRAY = "#3b6ea8", "#c0392b", "#7f7f7f"
LOG_METRICS = {"total_counts": "UMI counts per nucleus", "n_genes_by_counts": "genes detected per nucleus"}


def log10p(x):
    return np.log10(x + 1)


def inv_log10p(y):
    return np.power(10, y) - 1


def mad_bounds(values, n_mads):
    """Raw-MAD bounds at +/- n_mads from the median, as in qc.py (scipy's scale=1.0)."""
    med = np.nanmedian(values)
    mad = median_abs_deviation(values, nan_policy="omit")
    return med - n_mads * mad, med + n_mads * mad


def bounds_in_log(path):
    """The bounds qc.py logged, e.g. '3-MAD total_counts (log10): bounds=[266.7, 26119.5], flag 52'."""
    found = {}
    for line in open(path):
        match = re.search(r"MAD (\w+) \(log10\): bounds=\[([\d.]+), ([\d.]+)\]", line)
        if match:
            found[match.group(1)] = (float(match.group(2)), float(match.group(3)))
    return found


os.makedirs(args.outdir, exist_ok=True)
with open(args.log, "w") as log:
    def say(msg):
        print(f"[qc_report:{args.sample}] {msg}", file=log, flush=True)

    raw = sc.read_h5ad(args.raw)
    qcd = sc.read_h5ad(args.qc)
    say(f"raw: {raw.n_obs} cells; QC'd: {qcd.n_obs} cells")

    # The same metrics qc.py computes (it also flags ribo/hb genes, which do not change these three).
    raw.var["mt"] = raw.var_names.str.startswith(("mt-", "MT-"))
    sc.pp.calculate_qc_metrics(raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    above_floor = (raw.obs["total_counts"] >= args.min_counts).to_numpy()
    say(f"min {args.min_counts:g} counts: {int((~above_floor).sum())} below, {int(above_floor.sum())} above")

    # Bounds on the floored nuclei, as qc.py applies them; checked against the QC log when given.
    bounds = {}
    for metric in LOG_METRICS:
        lo, hi = mad_bounds(log10p(raw.obs[metric].to_numpy()[above_floor]), args.n_mads)
        bounds[metric] = (inv_log10p(lo), inv_log10p(hi))
        say(f"{args.n_mads:g}-MAD {metric}: bounds=[{bounds[metric][0]:.1f}, {bounds[metric][1]:.1f}]")
    if args.qc_log:
        logged = bounds_in_log(args.qc_log)
        for metric, (lo, hi) in bounds.items():
            if metric not in logged:
                sys.exit(f"{args.qc_log} has no bounds line for {metric}")
            if abs(logged[metric][0] - lo) > 0.06 or abs(logged[metric][1] - hi) > 0.06:
                sys.exit(f"{metric}: recomputed bounds [{lo:.1f}, {hi:.1f}] differ from the QC log's {logged[metric]}")
        say("bounds match the QC log")

    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    fig.suptitle(f"{args.sample} ({args.species}): {raw.n_obs:,} CellBender nuclei, {qcd.n_obs:,} after QC", y=0.995)

    for ax, (metric, label) in zip(axes[0], LOG_METRICS.items()):
        values = raw.obs[metric].to_numpy()
        positive = values[values > 0]
        bins = np.logspace(np.log10(positive.min()), np.log10(positive.max()), 60)
        ax.hist(positive, bins=bins, color=GRAY, alpha=0.8, label="all CellBender nuclei")
        ax.hist(values[above_floor], bins=bins, color=BLUE, alpha=0.8, label=f"above {args.min_counts:g}-count floor")
        lo, hi = bounds[metric]
        ax.axvline(lo, color=RED, linestyle="--", linewidth=1)
        ax.axvline(hi, color=RED, linestyle="--", linewidth=1, label=f"{args.n_mads:g}-MAD bounds [{lo:.0f}, {hi:.0f}]")
        if metric == "total_counts":
            ax.axvline(args.min_counts, color="black", linestyle=":", linewidth=1, label=f"floor {args.min_counts:g}")
        ax.set_xscale("log")
        ax.set_xlabel(label)
        ax.set_ylabel("nuclei")
        ax.legend(fontsize=8, frameon=False)

    ax = axes[1, 0]
    mito = raw.obs["pct_counts_mt"].to_numpy()[above_floor]
    ax.hist(mito, bins=np.linspace(0, max(mito.max(), args.max_pct_mito * 1.25), 60), color=BLUE, alpha=0.8)
    ax.axvline(args.max_pct_mito, color=RED, linestyle="--", linewidth=1,
               label=f"> {args.max_pct_mito:g}% removed ({int((mito > args.max_pct_mito).sum())})")
    ax.set_xlabel("mitochondrial counts (%), nuclei above the floor")
    ax.set_ylabel("nuclei")
    ax.set_yscale("log")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1, 1]
    if "doubletdetection_doublet_score" in qcd.obs:
        score = qcd.obs["doubletdetection_doublet_score"].to_numpy()
        ax.hist(score[np.isfinite(score)], bins=60, color=BLUE, alpha=0.8)
        ax.set_xlabel("DoubletDetection score, nuclei that passed QC")
    else:
        ax.text(0.5, 0.5, "no doublet scores in the QC'd object", ha="center", va="center", transform=ax.transAxes)
    ax.set_ylabel("nuclei")

    fig.tight_layout()
    figure_path = os.path.join(args.outdir, f"{args.sample}_qc.png")
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    say(f"wrote {figure_path}")

    with open(os.path.join(args.outdir, "summary.txt"), "w") as f:
        f.write(f"sample\t{args.sample}\nspecies\t{args.species}\n")
        f.write(f"cellbender_nuclei\t{raw.n_obs}\nabove_floor\t{int(above_floor.sum())}\n")
        for metric, (lo, hi) in bounds.items():
            f.write(f"{metric}_bounds\t{lo:.1f}\t{hi:.1f}\n")
        f.write(f"after_qc\t{qcd.n_obs}\npct_retained\t{100 * qcd.n_obs / raw.n_obs:.1f}\n")
    say("done")
