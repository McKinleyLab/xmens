#!/usr/bin/env python
"""Place heatmap-box SVGs side by side (left to right) with a shared colorbar SVG on the right,
optionally with a title (and a second subtitle line) above each box and gene names down the left side.

Each input is copied unchanged into its own <g id="..." transform="translate(x, y)"> group, so every
box keeps its exact size and cell colours; inner ids are prefixed with the panel id. Units are mm.
The colorbar bar is vertically centred on the boxes. Labels are live text (Helvetica Neue); their
extents are measured with matplotlib's font metrics to size the canvas.

usage: compose_heatmaps_svg.py OUT.svg --boxes ID=FILE [...] --colorbar FILE [--gap 2] [--cbar-gap 2]
           [--titles T1 T2 ...] [--subtitles S1 S2 ...] [--genes GENE_ORDER_TSV]
           [--title-pt 6] [--sub-pt 5] [--gene-pt 5] [--label-gap 0.8]
"""
import argparse
import copy
import xml.etree.ElementTree as ET

SVG = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG)
ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--boxes", nargs="+", required=True, help="id=path, left to right")
ap.add_argument("--colorbar", required=True)
ap.add_argument("--gap", type=float, default=2.0, help="mm between boxes")
ap.add_argument("--cbar-gap", type=float, default=2.0, help="mm between the last box and the colorbar")
ap.add_argument("--titles", nargs="*", default=None)
ap.add_argument("--subtitles", nargs="*", default=None)
ap.add_argument("--genes", default=None, help="TSV with a 'gene' column, in row order (top to bottom)")
ap.add_argument("--title-pt", type=float, default=6.0)
ap.add_argument("--sub-pt", type=float, default=5.0)
ap.add_argument("--gene-pt", type=float, default=5.0)
ap.add_argument("--label-gap", type=float, default=0.8, help="mm between labels and the boxes")
a = ap.parse_args()

PT = 25.4 / 72
FF = "'Helvetica Neue', Helvetica, Arial, sans-serif"
fmt = lambda v: f"{v:.5f}".rstrip("0").rstrip(".")
size = lambda root: tuple(float(root.get(k).removesuffix("mm")) for k in ("width", "height"))
if a.titles or a.subtitles or a.genes:
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    def ext(s, pt, italic=False):
        return TextPath((0, 0), s, size=pt * PT, prop=FontProperties(family="Helvetica Neue", style="italic" if italic else "normal")).get_extents()

boxes = []
for spec in a.boxes:
    gid, path = spec.split("=", 1)
    root = ET.parse(path).getroot()
    boxes.append((gid, root, *size(root)))
H = max(h for *_, h in boxes)
assert all(abs(h - H) < 1e-9 for *_, h in boxes), "boxes differ in height"
for lab, vals in (("titles", a.titles), ("subtitles", a.subtitles)):
    assert vals is None or len(vals) == len(boxes), f"need one of --{lab} per box"

cb = ET.parse(a.colorbar).getroot(); cbw, cbh = size(cb)
bar = list(cb.iter(f"{{{SVG}}}rect"))
bar_y0 = min(float(r.get("y")) for r in bar)
bar_h = max(float(r.get("y")) + float(r.get("height")) for r in bar) - bar_y0

x = 0.0; placed = []
for gid, root, w, h in boxes:
    placed.append((gid, root, x, 0.0, w)); x += w + a.gap
x += a.cbar_gap - a.gap
cb_y = (H - bar_h) / 2 - bar_y0
placed.append(("colorbar", cb, x, cb_y, cbw))
right = x + cbw

# ---- label geometry ----------------------------------------------------------------------
labels = []                                       # (id, attrs, text)
left, top = 0.0, min(0.0, cb_y)
if a.genes:
    import pandas as pd
    genes = list(pd.read_csv(a.genes, sep="\t").gene)
    nr = len(genes); pitch = H / nr
    cap = ext("H", a.gene_pt).height
    widest = max(ext(g, a.gene_pt, italic=True).x1 for g in genes)
    left = min(left, -(a.label_gap + widest + 0.3))
    for i, g in enumerate(genes):
        labels.append(("genes", {"x": fmt(-a.label_gap), "y": fmt((i + 0.5) * pitch + cap / 2), "text-anchor": "end",
                                 "font-style": "italic", "font-size": fmt(a.gene_pt * PT)}, g))
if a.titles or a.subtitles:
    y = -a.label_gap                              # baselines stack upward from just above the boxes
    for lab, vals, pt in (("subtitles", a.subtitles, a.sub_pt), ("titles", a.titles, a.title_pt)):
        if not vals:
            continue
        desc = max(-ext(v, pt).y0 for v in vals)  # room for descenders below the baseline
        y -= desc
        for (gid, _, gx, _, w), v in zip(placed, vals):
            assert ext(v, pt).width <= w + a.gap - 0.2, f"'{v}' is wider than its panel at {pt} pt"
            labels.append((lab, {"x": fmt(gx + w / 2), "y": fmt(y), "text-anchor": "middle", "font-size": fmt(pt * PT)}, v))
        y -= max(ext(v, pt).y1 for v in vals) + 0.4   # ascent + line gap
    top = min(top, y - 0.2)
bottom = max(H, cb_y + cbh)

out = ET.Element(f"{{{SVG}}}svg", {"width": f"{fmt(right - left)}mm", "height": f"{fmt(bottom - top)}mm",
                                   "viewBox": f"{fmt(left)} {fmt(top)} {fmt(right - left)} {fmt(bottom - top)}"})
out.append(ET.Comment(f" panels left to right: {', '.join(p[0] for p in placed)}; gap {a.gap:g} mm; units mm; "
                      f"boxes start at x = 0, y = 0 "))
for gid, root, gx, gy, _ in placed:
    g = ET.SubElement(out, f"{{{SVG}}}g", {"id": gid, "transform": f"translate({fmt(gx)},{fmt(gy)})"})
    for child in root:
        c = copy.deepcopy(child)
        for el in c.iter():                       # keep ids unique across panels
            if el.get("id"): el.set("id", f"{gid}_{el.get('id')}")
        g.append(c)
for grp in dict.fromkeys(l[0] for l in labels):
    g = ET.SubElement(out, f"{{{SVG}}}g", {"id": grp, "font-family": FF, "fill": "#000000"})
    for _, attrs, text in (l for l in labels if l[0] == grp):
        ET.SubElement(g, f"{{{SVG}}}text", attrs).text = text
ET.ElementTree(out).write(a.out, encoding="UTF-8", xml_declaration=True)
print({"out": a.out, "canvas_mm": (round(right - left, 3), round(bottom - top, 3)), "viewBox_origin": (round(left, 3), round(top, 3)),
       "box_x_mm": [round(p[2], 3) for p in placed[:-1]], "colorbar_xy_mm": (round(x, 3), round(cb_y, 3)), "n_labels": len(labels)})
