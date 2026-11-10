#!/usr/bin/env python
"""Render the 3 dpi pseudotime heatmap as a bare box (no outline, ticks or labels) plus a
matching standalone colorbar, both as SVG.

Heatmap: every cell is a vector <rect>; the SVG canvas *is* the heatmap, exactly BOX_W x BOX_H mm.
Colorbar: the same 256-level colormap lookup as the heatmap cells, drawn as vector rects, with
integer ticks and live text (Helvetica Neue) so labels stay editable in Illustrator.

Heatmap recipe (METHODS_pseudotime_heatmap.md): cells (--cells: one timepoint, 'induced', or 'all') ordered by the global no-apoptotic
RTK pseudotime (rank-normalised), 100 equal-count bins, mean log1p CP10K per bin,
gaussian_filter1d(sigma=2, axis=1, mode='reflect'), rows in published order, z-scored per gene.

usage: render_heatmap_box_svg.py PERCELL_TSV_GZ GENE_ORDER_TSV --box W H --zlim L [--outline-pt 0]
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties, findfont
from scipy.ndimage import gaussian_filter1d

ap = argparse.ArgumentParser()
ap.add_argument("percell"); ap.add_argument("order")
ap.add_argument("--box", nargs=2, type=float, default=[30.0, 54.786], metavar=("W_MM", "H_MM"))
ap.add_argument("--zlim", type=float, default=4.0, help="symmetric colour limit for z")
ap.add_argument("--scale", choices=["z", "minmax"], default="z",
                help="z: diverging RdBu_r over ±zlim; minmax: values already 0–1, white-to-red (upper half of RdBu_r) over 0..1")
ap.add_argument("--outline-pt", type=float, default=0.0, help="black outline around the box; 0 = none")
ap.add_argument("--cbar-mm", nargs=2, type=float, default=None, metavar=("W_MM", "H_MM"),
                help="colorbar bar size; default 1.6 mm x half the box height (as in the full figure)")
ap.add_argument("--font-pt", type=float, default=3.6)
ap.add_argument("--label", default="z per gene")
ap.add_argument("--cells", default="3_dpi", help="a timepoint label, 'induced' (all but uninduced) or 'all'")
ap.add_argument("--prefix", default=None, help="output prefix; default heatmap_<cells>")
ap.add_argument("--matrix", default=None, help="precomputed genes x columns z matrix (TSV, first column = gene), e.g. the Monocle-smoothed output of monocle_heatmap_repro.R; replaces the binned + Gaussian computation")
a = ap.parse_args()

BOX_W, BOX_H = a.box
N_BINS, SIGMA = 100, 2.0
PT = 25.4 / 72                                  # mm per point
EPS = 0.01                                      # mm overlap into the next-drawn cell; hides viewer seams
CMAP = plt.get_cmap("RdBu_r")                   # 256-entry lookup, same as matplotlib imshow
if a.scale == "minmax":                         # same white-to-red map as the earlier min–max deliverables
    CMAP = mpl.colors.LinearSegmentedColormap.from_list("rdbu_reds", plt.get_cmap("RdBu_r")(np.linspace(0.5, 1.0, 256)))
fmt = lambda v: f"{v:.5f}".rstrip("0").rstrip(".")
tag = f"{BOX_W:g}x{BOX_H:g}mm".replace(".", "p")
ztag = f"z{a.zlim:g}".replace(".", "p") if a.scale == "z" else "minmax"

# ---- data -------------------------------------------------------------------------------
pc = pd.read_csv(a.percell, sep="\t")
order = list(pd.read_csv(a.order, sep="\t").gene)
d = pc if a.cells == "all" else pc[pc.timepoint != "uninduced"] if a.cells == "induced" else pc[pc.timepoint == a.cells]
assert len(d), f"no cells for --cells {a.cells}"
PREFIX = a.prefix or f"heatmap_{a.cells.replace('_', '')}"
if a.matrix:                                   # precomputed matrix: rows reordered to the published gene order
    Zdf = pd.read_csv(a.matrix, sep="\t", index_col=0)
    assert set(order) <= set(Zdf.index), f"matrix lacks genes: {sorted(set(order) - set(Zdf.index))}"
    Z = Zdf.loc[order].to_numpy(float)
else:
    ix = np.argsort(d.pt_rank_norm.to_numpy(), kind="stable")
    X = d[[f"log_{g}" for g in order]].to_numpy()
    A = gaussian_filter1d(np.vstack([X[b].mean(axis=0) for b in np.array_split(ix, N_BINS)]).T, sigma=SIGMA, axis=1)
    Z = (A - A.mean(axis=1, keepdims=True)) / (A.std(axis=1, keepdims=True) + 1e-12)
if a.scale == "z":
    assert np.abs(Z).max() <= a.zlim, f"max |z| {np.abs(Z).max():.3f} exceeds ±{a.zlim}; cells would clip"
    norm = mpl.colors.Normalize(-a.zlim, a.zlim)
else:
    assert Z.min() >= -1e-12 and Z.max() <= 1 + 1e-12, "minmax matrix must lie in [0, 1]"
    norm = mpl.colors.Normalize(0.0, 1.0)


# ---- heatmap box --------------------------------------------------------------------------
def write_box(path):
    nr, nc = Z.shape
    cw, ch = BOX_W / nc, BOX_H / nr
    sw = a.outline_pt * PT
    pad = sw / 2 + 0.05 if sw > 0 else 0.0
    W, H = BOX_W + 2 * pad, BOX_H + 2 * pad
    rgba = CMAP(norm(Z))
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           f'<svg xmlns="http://www.w3.org/2000/svg" width="{fmt(W)}mm" height="{fmt(H)}mm" viewBox="0 0 {fmt(W)} {fmt(H)}">',
           f"<!-- heatmap {BOX_W:g} x {BOX_H:g} mm; {nr} genes x {nc} pseudotime bins; "
           + (f"z per gene, colour limits ±{a.zlim:g}" if a.scale == "z" else "min–max scaled per gene, 0 to 1") + " -->",
           '<g id="heatmap" shape-rendering="crispEdges" stroke="none">']
    for r in range(nr):
        for c in range(nc):
            w = cw + (EPS if c < nc - 1 else 0.0)
            h = ch + (EPS if r < nr - 1 else 0.0)
            out.append(f'<rect x="{fmt(pad + c * cw)}" y="{fmt(pad + r * ch)}" width="{fmt(w)}" height="{fmt(h)}" '
                       f'fill="{mpl.colors.to_hex(rgba[r, c])}"/>')
    out.append("</g>")
    if sw > 0:
        out.append(f'<rect id="outline" x="{fmt(pad)}" y="{fmt(pad)}" width="{fmt(BOX_W)}" height="{fmt(BOX_H)}" '
                   f'fill="none" stroke="#000000" stroke-width="{fmt(sw)}"/>')
    out.append("</svg>")
    open(path, "w").write("\n".join(out) + "\n")


# ---- colorbar -----------------------------------------------------------------------------
def write_cbar(path):
    bw, bh = a.cbar_mm if a.cbar_mm else (1.6, BOX_H * 0.5)
    fs = a.font_pt * PT                                         # font size in mm
    tick_len, tick_w, tick_pad, lab_pad = 0.6 * PT, 0.25 * PT, 0.8 * PT, 1.2 * PT
    prop = FontProperties(family="Helvetica Neue")
    font_file = findfont(prop)
    ext = lambda s: TextPath((0, 0), s, size=fs, prop=prop).get_extents()
    if a.scale == "z":
        ticks = np.arange(np.ceil(-a.zlim), np.floor(a.zlim) + 1).astype(int)
        labels = [("−" if t < 0 else "") + str(abs(t)) for t in ticks]
    else:
        ticks, labels = [0.0, 0.5, 1.0], ["0", "0.5", "1"]
    cap = ext("0").height                                       # digit height, for vertical centring
    pad = cap / 2 + 0.3                                         # room for end labels above/below the bar
    x0, y0 = 0.05, pad
    tx = x0 + bw + tick_len + tick_pad                          # left edge of tick labels
    lab_w = max(ext(s).x1 for s in labels)
    lx = tx + lab_w + lab_pad + ext(a.label).y1                 # baseline x of the rotated label (ascent to its left)
    W, H = lx + abs(ext(a.label).y0) + 0.3, bh + 2 * pad
    yv = lambda v: y0 + bh * (1 - norm(v))                      # value -> y (top = +zlim)
    n = CMAP.N
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           f'<svg xmlns="http://www.w3.org/2000/svg" width="{fmt(W)}mm" height="{fmt(H)}mm" viewBox="0 0 {fmt(W)} {fmt(H)}">',
           "<!-- colorbar for " + (f"z per gene, ±{a.zlim:g}" if a.scale == "z" else "min–max scaled per gene, 0 to 1")
           + f"; bar {bw:g} x {bh:g} mm; {n} colour levels -->",
           '<g id="colorbar" shape-rendering="crispEdges" stroke="none">']
    for i in range(n):                                          # level i spans values [i, i+1)/n of the range
        top = y0 + bh * (1 - (i + 1) / n)
        h = bh / n + (EPS if i > 0 else 0.0)                    # overlap downward into the previously drawn level
        out.append(f'<rect x="{fmt(x0)}" y="{fmt(top)}" width="{fmt(bw)}" height="{fmt(h)}" '
                   f'fill="{mpl.colors.to_hex(CMAP(i))}"/>')
    out.append("</g>")
    out.append(f'<g id="ticks" stroke="#666666" stroke-width="{fmt(tick_w)}">')
    for t in ticks:
        out.append(f'<line x1="{fmt(x0 + bw)}" y1="{fmt(yv(t))}" x2="{fmt(x0 + bw + tick_len)}" y2="{fmt(yv(t))}"/>')
    out.append("</g>")
    ff = "font-family=\"'Helvetica Neue', Helvetica, Arial, sans-serif\""
    out.append(f'<g id="tick-labels" {ff} font-size="{fmt(fs)}" fill="#000000">')
    for t, s in zip(ticks, labels):
        out.append(f'<text x="{fmt(tx)}" y="{fmt(yv(t) + cap / 2)}">{s}</text>')
    out.append("</g>")
    out.append(f'<text id="label" {ff} font-size="{fmt(fs)}" fill="#000000" text-anchor="middle" '
               f'transform="translate({fmt(lx)},{fmt(y0 + bh / 2)}) rotate(-90)">{a.label}</text>')
    out.append("</svg>")
    open(path, "w").write("\n".join(out) + "\n")
    return dict(bar_mm=(bw, round(bh, 4)), canvas_mm=(round(W, 3), round(H, 3)), ticks=list(ticks),
                font=font_file.split("/")[-1], bar_x0=x0, bar_y0=y0)


box_path = f"{PREFIX}_{ztag}_box_{tag}.svg"
cb_path = f"{PREFIX}_{ztag}_colorbar_{tag}.svg"
write_box(box_path)
info = write_cbar(cb_path)
print(box_path, cb_path, info, sep="\n")
print(f"cells={len(d)} z_range=({Z.min():.3f},{Z.max():.3f}) "
      + (f"red_range_used={100 * Z.max() / a.zlim:.1f}% blue_range_used={100 * -Z.min() / a.zlim:.1f}%" if a.scale == "z" else "scale=minmax 0..1"))
