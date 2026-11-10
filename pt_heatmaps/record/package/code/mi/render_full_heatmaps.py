#!/usr/bin/env python
"""Full pseudotime heatmaps: every top-decile pseudotime-associated gene, one panel per species, rows sorted by
peak pseudotime. Remake of the first-submission Fig. 4C with the current pseudotime, MI arms and Monocle smoother.

Gene sets (Methods, "Pseudotime-associated genes and cross-species conservation"): genes with MI above the 90th
percentile of MI among each species' significant genes, restricted to the background of strict one-to-one orthologs
significant in both species (final objects + MI of record: 1,201 human, 1,133 mouse; 373 shared).
Curves: monocle_conserved_smooth.R (Monocle 2 plot_pseudotime_heatmap smoother: negative binomial, ns(df = 3),
100 evenly spaced points of rank-normalized pseudotime). Transform as in plot_pseudotime_heatmap: log10(curve + 1),
z per gene (sd, n - 1), clipped at +/- zclip.
Row order: index of the curve maximum (peak pseudotime). Ties are the genes peaking at either end of pseudotime;
genes peaking at the start are ordered by when their z first falls below 0, genes peaking at the end by when it
last rises above 0, so each end block runs from the sharpest to the most gradual transition. (A centre-of-mass
tie-break stacked the few genes with a second hump at the block edges, where they showed as horizontal seams.)

usage: render_full_heatmaps.py --curves-human TSV --curves-mouse TSV --genes-human CSV --genes-mouse CSV
                               --cells-human TSV --cells-mouse TSV --out PREFIX [--zclip 3]
(the cells tables are only used for the nucleus counts in the panel titles)
Outputs: PREFIX.svg/.pdf/.png, PREFIX_gene_order_{human,mouse}.csv, PREFIX_provenance.json
"""
import argparse, datetime, hashlib, json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import ListedColormap, Normalize

ap = argparse.ArgumentParser()
for name in ("--curves-human", "--curves-mouse", "--genes-human", "--genes-mouse", "--cells-human", "--cells-mouse", "--out"):
    ap.add_argument(name, required=True)
ap.add_argument("--zclip", type=float, default=3.0)
ap.add_argument("--shared-column", action="store_true",
                help="draw a black/white column left of each panel marking genes also in the other species' top set "
                     "(off by default; the shared flag is always written to the gene-order tables)")
ap.add_argument("--order", choices=("peak", "centroid"), default="peak",
                help="peak: curve maximum, end-of-pseudotime ties by crossing time; centroid: z-weighted centre of mass "
                     "(the ordering of the 28-gene cross-species heatmap)")
args = ap.parse_args()

HELVETICA = "/System/Library/Fonts/HelveticaNeue.ttc"
if os.path.exists(HELVETICA):
    fm.fontManager.addfont(HELVETICA)
# byte-reproducible vector output: fixed SVG element ids, no creation dates
plt.rcParams["svg.hashsalt"] = "xmens-pseudotime-heatmaps"
STABLE_META = {"svg": {"Date": None}, "pdf": {"CreationDate": None}}
plt.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "DejaVu Sans"],
                     "svg.fonttype": "none", "pdf.fonttype": 42, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "xtick.major.size": 2.0, "xtick.major.pad": 1.5})
SIZE_BASE, SIZE_SMALL, SIZE_TICK = 8, 7, 6

def ordered_z(curves_tsv, genes_csv):
    """Monocle z-scores for every gene in the set, rows sorted by peak pseudotime then centroid."""
    curves = pd.read_csv(curves_tsv, sep="\t", index_col=0)
    fit_used = curves.pop("fit")
    grid = pd.read_csv(curves_tsv.replace(".tsv", "_grid.tsv"), sep="\t")["Pseudotime"].to_numpy()
    genes = pd.read_csv(genes_csv)
    missing = set(genes.gene) - set(curves.index)
    assert not missing, f"{curves_tsv} lacks {len(missing)} genes"
    log_curve = np.log10(curves.loc[genes.gene].to_numpy(float) + 1)
    z = (log_curve - log_curve.mean(1, keepdims=True)) / log_curve.std(1, ddof=1, keepdims=True)
    peak = z.argmax(1)
    below = z < 0
    # first column where z < 0 (start-peaking genes) and first column after the last z < 0 (end-peaking genes)
    fall = np.where(below.any(1), below.argmax(1), z.shape[1])
    rise = np.where(below.any(1), z.shape[1] - below[:, ::-1].argmax(1), 0)
    positive = np.clip(z, 0, None); t = (grid - grid.min()) / (grid.max() - grid.min())
    centroid = (positive * t).sum(1) / positive.sum(1)
    tie_break = np.where(peak == 0, fall, np.where(peak == z.shape[1] - 1, rise, centroid))
    table = genes.assign(peak_index=peak, peak_pseudotime=grid[peak], fall_index=fall, rise_index=rise,
                         centroid=centroid, fit=fit_used.reindex(genes.gene).to_numpy())
    order = np.lexsort((tie_break, peak)) if args.order == "peak" else np.lexsort((peak, centroid))
    return z[order], table.iloc[order].reset_index(drop=True), int((np.abs(z) > args.zclip).sum())


panels = {}
for species in ("human", "mouse"):
    z, table, n_clipped = ordered_z(getattr(args, f"curves_{species}"), getattr(args, f"genes_{species}"))
    n_cells = len(pd.read_csv(getattr(args, f"cells_{species}"), sep="\t", usecols=["cell"]))
    panels[species] = dict(z=z, table=table, n_clipped=n_clipped, n_cells=n_cells)

# ---- layout (mm from the top-left corner) ------------------------------------------------------------------------
MM = 1 / 25.4
W_MM, H_MM = 160.0, 156.0
fig = plt.figure(figsize=(W_MM * MM, H_MM * MM))
def axes_mm(x, y_top, w, h):
    return fig.add_axes([x / W_MM, 1 - (y_top + h) / H_MM, w / W_MM, h / H_MM])

X0 = {"human": 16.0, "mouse": 80.0}
HEAT_W, HEAT_TOP, HEAT_H = 52.0, 8.0, 136.0
SHARED_W, SHARED_GAP = 1.6, 0.6
norm = Normalize(-args.zclip, args.zclip); cmap = plt.get_cmap("RdBu_r")

for species, p in panels.items():
    x0, n_genes = X0[species], len(p["table"])
    heat = axes_mm(x0, HEAT_TOP, HEAT_W, HEAT_H)
    heat.imshow(np.clip(p["z"], -args.zclip, args.zclip), cmap=cmap, norm=norm, aspect="auto",
                interpolation="nearest", extent=(0, 1, n_genes, 0), rasterized=True)
    heat.set_xlim(0, 1); heat.set_ylim(n_genes, 0)
    heat.set_xticks([0, 0.5, 1]); heat.set_xticklabels(["0", "0.5", "1"], fontsize=SIZE_TICK)
    heat.set_yticks([])
    heat.set_xlabel("pseudotime (rank-normalized)", fontsize=SIZE_BASE, labelpad=2)
    heat.set_title(f"{species.capitalize()} ({p['n_cells']:,} nuclei)", loc="left", fontsize=SIZE_BASE, pad=3)

    order_text = "ordered by peak pseudotime" if args.order == "peak" else "ordered by pseudotime centroid"
    if not args.shared_column:
        heat.set_ylabel(f"{n_genes:,} genes, {order_text}", fontsize=SIZE_SMALL, labelpad=2)
        continue
    # optional: column marking the genes also in the other species' top set (the 1:1-ortholog overlap)
    shared = axes_mm(x0 - SHARED_GAP - SHARED_W, HEAT_TOP, SHARED_W, HEAT_H)
    shared.imshow(p["table"]["shared"].to_numpy(int)[:, None], cmap=ListedColormap(["white", "black"]), vmin=0, vmax=1,
                  aspect="auto", interpolation="nearest", extent=(0, 1, n_genes, 0), rasterized=True)
    shared.set_xticks([]); shared.set_yticks([])
    for spine in shared.spines.values():
        spine.set_linewidth(0.5)
    other = "mouse" if species == "human" else "human"
    shared.set_ylabel(f"{n_genes:,} genes, {order_text}\n(black: shared with {other}, "
                      f"{int(p['table']['shared'].sum())})", fontsize=SIZE_SMALL, labelpad=2)

colorbar_ax = axes_mm(138.0, 110.0, 3.0, 30.0)
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=colorbar_ax)
ticks = np.arange(-int(args.zclip), int(args.zclip) + 1)
cb.set_ticks(ticks); cb.set_ticklabels([("−" if t < 0 else "") + str(abs(t)) for t in ticks])
cb.ax.tick_params(labelsize=SIZE_TICK, width=0.5, length=2, pad=1.5); cb.outline.set_linewidth(0.5)
cb.set_label(f"expression, z per gene\n(clipped at ±{args.zclip:g})", fontsize=SIZE_SMALL)

for ext in ("svg", "pdf"):
    fig.savefig(f"{args.out}.{ext}", dpi=600, metadata=STABLE_META[ext])
fig.savefig(f"{args.out}.png", dpi=600)

# render check: no visible text boxes overlap each other or leave the canvas
fig.set_dpi(100); renderer = fig.canvas.get_renderer()
boxes = [(t.get_text(), t.get_window_extent(renderer)) for t in fig.findobj(matplotlib.text.Text)
         if t.get_text().strip() and t.get_visible()]
overlaps = [(a, b) for i, (a, ba) in enumerate(boxes) for b, bb in boxes[i + 1:] if ba.overlaps(bb)]
outside = [a for a, b in boxes if not fig.bbox.contains(b.x0, b.y0) or not fig.bbox.contains(b.x1, b.y1)]
plt.close(fig)

for species, p in panels.items():
    p["table"].assign(row=np.arange(1, len(p["table"]) + 1)).to_csv(f"{args.out}_gene_order_{species}.csv", index=False)
sha = lambda path: hashlib.sha256(open(path, "rb").read()).hexdigest()
inputs = {k: getattr(args, k) for k in ("curves_human", "curves_mouse", "genes_human", "genes_mouse", "cells_human", "cells_mouse")}
provenance = dict(written=datetime.datetime.now().isoformat(timespec="seconds"),
                  inputs={k: dict(path=os.path.abspath(v), sha256=sha(v)) for k, v in inputs.items()},
                  n_genes={s: len(p["table"]) for s, p in panels.items()},
                  n_shared={s: int(p["table"]["shared"].sum()) for s, p in panels.items()},
                  n_cells={s: p["n_cells"] for s, p in panels.items()},
                  clipped_values={s: p["n_clipped"] for s, p in panels.items()}, zclip=args.zclip,
                  shared_column=bool(args.shared_column), order=args.order,
                  fallback_fits={s: int((p["table"]["fit"] != "trend").sum()) for s, p in panels.items()})
json.dump(provenance, open(f"{args.out}_provenance.json", "w"), indent=1)
print(json.dumps(dict(n_genes=provenance["n_genes"], n_shared=provenance["n_shared"], clipped=provenance["clipped_values"], fallback=provenance["fallback_fits"], overlaps=overlaps, outside=outside), default=str))
