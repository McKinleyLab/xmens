"""Compare a run of this pipeline with the heatmaps of record, step by step, for each step the run has finished.

  1 inputs   the nuclei, pseudotime, labels and timepoints, and the gene lists (MI top genes, Fig. 5F genes, Yang
             panel), against the inputs of the heatmap package of record; the counts of the Fig. 4H and Fig. 5F genes,
             nucleus by nucleus, against the package's per-nucleus count tables
  2 monocle  Monocle's dispersion table (genes, mean, empirical dispersion) and dispersion trend against the package's
             R port, which used the sample variance (n - 1) where Monocle uses the population variance (n): the port's
             dispersions are converted to the population variance (with the size factors from step 1's counts) and must
             then match Monocle's; fallback fits; Monocle's z matrices, clipped at +/-3, against the package's z of
             Fig. 4H and Fig. 5F
  3 figures  the MI-heatmap gene order against the package's; each SVG against the figure of record: byte for byte,
             else its heatmap cell colours, embedded images, text strings and the rest of the file (positions and
             layout) separately
Writes one row per comparison to --out (TSV) and prints a summary.
"""
import argparse
import base64
import difflib
import hashlib
import io
import re
from pathlib import Path

import matplotlib.image
import numpy as np
import pandas as pd
import scipy.io
from scipy.stats import spearmanr

SPECIES = ("mouse", "human")
FIGURES = [  # (what, run file under figures/, record file under record/package/)
    ("Fig. 4H box", "fig4h/heatmap_monocle_induced_z3_box_30x59p887mm.svg",
     "yang_induced/heatmap_monocle_induced_z3_box_30x59p887mm.svg"),
    ("Fig. 4H colorbar, 5 pt", "fig4h/heatmap_monocle_z3_colorbar_5pt_30x59p887mm.svg",
     "yang_induced/heatmap_monocle_z3_colorbar_5pt_30x59p887mm.svg"),
    ("Fig. 4H labelled", "fig4h/heatmap_monocle_induced_labelled_30x59p887mm.svg",
     "yang_induced/heatmap_monocle_induced_labelled_30x59p887mm.svg"),
    ("Fig. 5F", "fig5f/crossspecies28_monocle_clip3_final_objects_nostrips.svg",
     "crossspecies/crossspecies28_monocle_clip3_final_objects_nostrips.svg"),
    ("MI heatmaps", "mi/mi_top_genes_heatmap_centroid_order.svg",
     "mi_heatmaps/mi_top_genes_heatmap_centroid_order.svg"),
]
Z_CLOSE = 0.01  # largest |difference| in z (colour limits +/-3) reported as "close"

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--run", type=Path, required=True, help="the run folder ($OUT)")
parser.add_argument("--record", type=Path, required=True, help="the record/ folder")
parser.add_argument("--out", type=Path, required=True, help="output TSV; must not exist yet")
args = parser.parse_args()
if args.out.exists():
    raise SystemExit(f"{args.out} exists; not overwriting")
package = args.record / "package"

rows = []


def report(step, species, what, result, detail=""):
    rows.append({"step": step, "species": species, "what": what, "result": result, "detail": detail})
    print(f"{step} {species} {what}: {result}{'; ' + detail if detail else ''}", flush=True)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def same_text_rows(run_file, record_file, **read_options):
    """'identical' if the files are equal byte for byte, or hold the same table as text, row for row (a gzip file also
    stores its time stamp); 'same rows, other order' if they hold the same rows as text in another order; otherwise
    'DIFFERS'. Also returns the run table as text."""
    if sha256(run_file) == sha256(record_file):
        return "identical", pd.read_csv(run_file, dtype=str, keep_default_na=False, **read_options)
    run = pd.read_csv(run_file, dtype=str, keep_default_na=False, **read_options)
    record = pd.read_csv(record_file, dtype=str, keep_default_na=False, **read_options)
    if run.equals(record):
        return "identical", run
    if list(run.columns) == list(record.columns) and \
            sorted(map(tuple, run.values.tolist())) == sorted(map(tuple, record.values.tolist())):
        return "same rows, other order", run
    return "DIFFERS", run


def compare_z(run_file, record_file, sep, genes=None):
    """Monocle's z (step 2) against a z matrix of the record (genes x pseudotime points), both clipped at +/-3 as when
    drawn (the record's Fig. 4H z was stored clipped, its Fig. 5F z unclipped): the same genes in the same order, and
    the largest difference."""
    run = pd.read_csv(run_file, sep="\t", index_col=0).drop(columns="fit").clip(-3, 3)
    record = pd.read_csv(record_file, sep=sep, index_col=0)
    record = record.loc[:, ~record.columns.str.startswith("Unnamed")].clip(-3, 3)   # a trailing tab adds a column
    if genes is not None:
        run = run.loc[genes]
    if list(run.index) != list(record.index) or run.shape != record.shape:
        same_genes = sorted(run.index) == sorted(record.index)
        return "DIFFERS", (f"genes {'in another order' if same_genes else 'differ'}: {len(run)} x {run.shape[1]}, "
                           f"record {len(record)} x {record.shape[1]}")
    difference = np.abs(run.to_numpy() - record.to_numpy())
    largest = difference.max()
    gene = run.index[int(difference.max(axis=1).argmax())]
    result = "identical" if largest == 0 else ("close" if largest <= Z_CLOSE else "DIFFERS")
    detail = (f"{len(run)} genes x {run.shape[1]} points, same genes in the same order; largest |difference| in z "
              f"{largest:.2g} ({gene}), mean {difference.mean():.2g}; clipped at +/-3: "
              f"{int((run.abs() >= 3).sum().sum())}"
              f", record {int((record.abs() >= 3).sum().sum())}")
    return result, detail


# ---- SVG parts -------------------------------------------------------------------------------------------------------
FILL = re.compile(r'(fill(?:="|:\s*))(#[0-9a-fA-F]{6})')
IMAGE = re.compile(r'(xlink:href="data:image/png;base64,)([^"]+)(")')
TEXT = re.compile(r"(<text[^>]*>)([^<]*)(</text>)")
MATPLOTLIB = re.compile(r"<dc:title>Matplotlib v[^<]*</dc:title>")


def svg_parts(path):
    svg = Path(path).read_text()
    fills = [m.group(2).lower() for m in FILL.finditer(svg)]
    images = [matplotlib.image.imread(io.BytesIO(base64.b64decode(re.sub(r"\s", "", m.group(2)))), format="png")
              for m in IMAGE.finditer(svg)]
    texts = [m.group(2) for m in TEXT.finditer(svg)]
    rest = MATPLOTLIB.sub("<dc:title>Matplotlib</dc:title>", svg)            # the version that drew it
    rest = TEXT.sub(r"\1<text>\3", IMAGE.sub(r"\1<png>\3", FILL.sub(r"\1#colour", rest)))
    return fills, images, texts, rest.splitlines()


COLOUR_CLOSE = 3  # largest change of a colour channel (of 255) reported as "close"


def compare_svg(run_file, record_file):
    """'identical' byte for byte; 'same drawing' if only the matplotlib version differs; 'close' if, beyond that, only
    colours differ, by at most COLOUR_CLOSE of 255 in a channel; 'same content, layout differs' if the colours (to
    COLOUR_CLOSE) and text strings are the same but positions or the size of an embedded image differ, as when the
    text is measured with another font; otherwise 'DIFFERS'."""
    if sha256(run_file) == sha256(record_file):
        return "identical", "byte for byte"
    run, record = svg_parts(run_file), svg_parts(record_file)
    parts = []
    colour_change = 0       # largest change of a colour channel in the heatmap cells and images
    other_change = False    # content: the number of cells or images, or the text strings
    layout_change = False   # positions, or the pixel size of an embedded image
    fills, record_fills = run[0], record[0]
    if len(fills) != len(record_fills):
        parts.append(f"colours: {len(fills):,} fills, record {len(record_fills):,}")
        other_change = True
    else:
        differ = [(a, b) for a, b in zip(fills, record_fills) if a != b]
        if differ:
            largest = max(max(abs(int(a[i:i + 2], 16) - int(b[i:i + 2], 16)) for i in (1, 3, 5)) for a, b in differ)
            colour_change = max(colour_change, largest)
            parts.append(f"colours: {len(differ):,} of {len(fills):,} fills differ, by at most {largest} of 255 in a "
                         f"channel")
        else:
            parts.append(f"colours: all {len(fills):,} fills identical")
    for i, (image, record_image) in enumerate(zip(run[1], record[1])):
        if image.shape != record_image.shape:
            parts.append(f"image {i + 1}: {image.shape[1]} x {image.shape[0]} px, record "
                         f"{record_image.shape[1]} x {record_image.shape[0]} (not compared)")
            layout_change = True
            continue
        differ = (image != record_image).any(axis=-1) if image.ndim == 3 else image != record_image
        if differ.any():
            largest = int(round(np.abs(image - record_image).max() * 255))
            colour_change = max(colour_change, largest)
            parts.append(f"image {i + 1}: {int(differ.sum()):,} of {differ.size:,} px differ, by at most {largest} "
                         f"of 255")
        else:
            parts.append(f"image {i + 1}: all {differ.size:,} px identical")
    if len(run[1]) != len(record[1]):
        parts.append(f"{len(run[1])} images, record {len(record[1])}")
        other_change = True
    if run[2] == record[2]:
        parts.append(f"text: the same {len(run[2])} strings")
    else:
        changed = [f"{a!r} (record {b!r})" for a, b in zip(run[2], record[2]) if a != b][:3]
        parts.append(f"text: {len(run[2])} strings, record {len(record[2])}; differing: {', '.join(changed)}")
        other_change = True
    differ = [line for line in difflib.unified_diff(record[3], run[3], lineterm="", n=0)
              if line[:1] in "+-" and not line.startswith(("+++", "---"))]
    if differ:
        parts.append(f"rest of the file (positions, layout): {sum(line.startswith('+') for line in differ):,} of "
                     f"{len(run[3]):,} lines differ")
        layout_change = True
    else:
        parts.append("rest of the file identical apart from the matplotlib version")
    if other_change or colour_change > COLOUR_CLOSE:
        result = "DIFFERS"
    elif layout_change:
        result = "same content, layout differs"
    else:
        result = "same drawing" if colour_change == 0 else "close"
    return result, "; ".join(parts)


# ---- 1 inputs --------------------------------------------------------------------------------------------------------
inputs = args.run / "inputs"
poisson_term = {}  # species -> (mean of 1 / size factor, number of nuclei), from the counts
if (inputs / "cells_mouse.tsv.gz").exists():
    for species in SPECIES:
        result, cells = same_text_rows(inputs / f"cells_{species}.tsv.gz",
                                       package / "inputs" / "common" / f"cells_{species}.tsv.gz", sep="\t")
        report("1 inputs", species, "nuclei, pseudotime, labels, timepoints", result,
               f"{len(cells):,} nuclei, compared as text (the gzip files store a time stamp)")
        result, top = same_text_rows(inputs / f"mi_top_{species}.csv",
                                     package / "inputs" / "mi" / f"topdecile_{species}.csv")
        report("1 inputs", species, "MI top genes", result, f"{len(top):,} genes")
    result, fig5f = same_text_rows(inputs / "fig5f_genes.csv",
                                   package / "inputs" / "crossspecies" / "conserved_profiles_selected_genes.csv")
    report("1 inputs", "both", "Fig. 5F genes", result, f"{len(fig5f)} genes")
    yang = pd.read_csv(inputs / "yang_panel.tsv", sep="\t")
    record_yang = pd.read_csv(package / "inputs" / "yang" / "gene_peak_order.tsv", sep="\t")
    same = yang["gene"].tolist() == record_yang["gene"].tolist()
    report("1 inputs", "mouse", "Yang et al. panel, published order", "identical" if same else "DIFFERS",
           f"{len(yang)} genes")

    # The counts of the drawn genes, nucleus by nucleus, against the package's per-nucleus tables.
    checks = {"mouse": [(yang["symbol"].tolist(), yang["gene"].tolist(), "yang/percell_genes_final.tsv.gz", "Fig. 4H"),
                        (fig5f["mouse"].tolist(), fig5f["mouse"].tolist(), "crossspecies/percell_mouse_sel28.tsv.gz",
                         "Fig. 5F")],
              "human": [(fig5f["human"].tolist(), fig5f["human"].tolist(), "crossspecies/percell_human_sel28.tsv.gz",
                         "Fig. 5F")]}
    for species in SPECIES:
        counts_file = inputs / species / "counts.mtx"
        if not counts_file.exists():
            continue
        counts = scipy.io.mmread(counts_file).tocsr()
        # Monocle's size factors (total counts over their geometric mean) give the Poisson term of its dispersion
        # estimate, the mean of 1 / size factor over the nuclei; step 2 uses it to convert the record's dispersions.
        total = np.asarray(counts.sum(axis=0)).ravel().astype(float)
        poisson_term[species] = (float(np.mean(np.exp(np.mean(np.log(total))) / total)), counts.shape[1])
        gene_row = {gene: i for i, gene in enumerate((inputs / species / "genes.txt").read_text().splitlines())}
        cell_names = pd.read_csv(inputs / f"cells_{species}.tsv.gz", sep="\t", usecols=["cell"])["cell"]
        report("1 inputs", species, "counts matrix", "read", f"{counts.shape[0]:,} genes x {counts.shape[1]:,} nuclei, "
               f"{counts.nnz:,} non-zero counts")
        for symbols, columns, record_file, figure in checks[species]:
            record = pd.read_csv(package / "inputs" / record_file, sep="\t",
                                 usecols=["cell"] + [f"cnt_{c}" for c in columns]).set_index("cell")
            ours = pd.DataFrame(counts[[gene_row[s] for s in symbols]].T.toarray(), index=cell_names,
                                columns=[f"cnt_{c}" for c in columns])
            if set(record.index) != set(ours.index):
                report("1 inputs", species, f"counts of the {figure} genes", "DIFFERS",
                       f"nuclei differ: {len(ours):,}, record {len(record):,}")
                continue
            n_differ = int((ours.loc[record.index].to_numpy() != record.to_numpy()).sum())
            report("1 inputs", species, f"counts of the {figure} genes", "identical" if n_differ == 0 else "DIFFERS",
                   f"{len(columns)} genes x {len(record):,} nuclei; {n_differ:,} counts differ")
        del counts

# ---- 2 monocle -------------------------------------------------------------------------------------------------------
z_dir = args.run / "z"
port_trend = pd.read_csv(args.record / "port_dispersion_trend.tsv", sep="\t").set_index("species")
for species in SPECIES:
    if not (z_dir / f"z_{species}.tsv").exists():
        continue
    monocle = pd.read_csv(z_dir / f"dispersion_{species}.tsv.gz", sep="\t").set_index("gene_id")
    record = pd.read_csv(package / "inputs" / "common" / f"disp_table_{species}.tsv.gz", sep="\t").set_index("gene")
    both = monocle.index.intersection(record.index)
    report("2 monocle", species, "dispersion table, genes", "identical" if len(both) == len(monocle) == len(record)
           else "DIFFERS", f"Monocle {len(monocle):,} genes, record {len(record):,}, in both {len(both):,}")
    mean = monocle.loc[both, "mean_expression"].to_numpy()
    record_mean = record.loc[both, "mu"].to_numpy()
    relative = np.abs(mean - record_mean) / record_mean
    report("2 monocle", species, "mean of size-factor-normalized counts",
           "identical" if relative.max() == 0 else ("close" if relative.max() < 1e-9 else "DIFFERS"),
           f"largest relative difference {relative.max():.1e}")
    # Monocle: (population variance - xim * mean) / mean^2, negatives set to 0; the record used the sample variance.
    disp = monocle.loc[both, "dispersion_empirical"].to_numpy()
    record_disp = record.loc[both, "disp"].to_numpy()
    positive = record_disp > 0
    relative = (disp[positive] - record_disp[positive]) / record_disp[positive]
    detail = (f"Monocle (population variance, n) lower than the record (sample variance, n - 1) for "
              f"{(relative < 0).mean():.1%} of the {positive.sum():,} genes with dispersion > 0, by a median "
              f"{-np.median(relative):.1e} of the value; dispersion 0 (negative estimate) in "
              f"{int((disp == 0).sum()):,} genes, record {int((record_disp == 0).sum()):,}")
    if species in poisson_term:
        xim, n = poisson_term[species]
        variance = record_disp * record_mean ** 2 + xim * record_mean        # the record's sample variance
        expected = np.where(positive, np.maximum(0, ((n - 1) / n * variance - xim * record_mean) / record_mean ** 2), 0)
        largest = np.abs(disp - expected).max()
        zero_mismatch = int(((disp == 0) != (expected == 0)).sum())
        result = "as expected" if largest < 1e-6 and zero_mismatch == 0 else "DIFFERS"
        detail += (f"; the record's values converted to the population variance (n = {n:,} nuclei, mean 1 / size "
                   f"factor {xim:.6g}) match Monocle's to {largest:.1e} (largest absolute difference), with "
                   f"{zero_mismatch} genes differing in being 0")
    else:
        result = "not converted"
        detail += "; not converted to the population variance: step 1's counts (for the size factors) are missing"
    report("2 monocle", species, "empirical dispersion", result, detail)
    trend = pd.read_csv(z_dir / f"dispersion_{species}_trend.tsv", sep="\t").set_index("coefficient")["value"]
    a, b = trend["asymptDisp"], trend["extraPois"]
    record_a, record_b = port_trend.loc[species, "asymptDisp"], port_trend.loc[species, "extraPois"]
    same4 = f"{a:.4g}" == f"{record_a:.4g}" and f"{b:.4g}" == f"{record_b:.4g}"
    report("2 monocle", species, "dispersion trend a + b / mean", "same to 4 digits" if same4 else "close"
           if max(abs(a / record_a - 1), abs(b / record_b - 1)) < 0.01 else "DIFFERS",
           f"Monocle a {a:.6g}, b {b:.6g}; record (printed to 4 digits) a {record_a:.4g}, b {record_b:.4g}")
    for name in [f"z_{species}"] + (["fig4h_z"] if species == "mouse" else []):
        fit = pd.read_csv(z_dir / f"{name}.tsv", sep="\t", usecols=["gene", "fit"])
        n_fallback = int((fit["fit"] != "trend").sum())
        report("2 monocle", species, f"{name}.tsv fits", "as record" if n_fallback == 0 else "DIFFERS",
               f"{len(fit):,} genes, {n_fallback} fallback fits (record: 0)")
    if species == "mouse":
        result, detail = compare_z(z_dir / "fig4h_z.tsv", package / "yang_induced" / "heatmap_monocle_induced_z.tsv",
                                   "\t")
        report("2 monocle", "mouse", "Fig. 4H z", result, detail)
    fig5f_genes = pd.read_csv(inputs / "fig5f_genes.csv")[species].tolist()
    result, detail = compare_z(z_dir / f"z_{species}.tsv", package / "crossspecies" /
                               f"crossspecies28_monocle_clip3_final_objects_nostrips_z_{species}.csv", ",", fig5f_genes)
    report("2 monocle", species, "Fig. 5F z", result, detail)

# ---- 3 figures -------------------------------------------------------------------------------------------------------
figures = args.run / "figures"
if figures.exists():
    for species in SPECIES:
        name = f"mi_top_genes_heatmap_centroid_order_gene_order_{species}.csv"
        run, record = pd.read_csv(figures / "mi" / name), pd.read_csv(package / "mi_heatmaps" / name)
        if sorted(run["gene"]) != sorted(record["gene"]):
            report("3 figures", species, "MI heatmap rows", "DIFFERS",
                   f"genes differ: {len(run)}, record {len(record)}")
            continue
        record_rows = run["gene"].map(dict(zip(record["gene"], range(len(record))))).to_numpy()
        move = np.abs(np.arange(len(run)) - record_rows)
        rho = spearmanr(np.arange(len(run)), record_rows).statistic
        merged = run.merge(record, on="gene", suffixes=("", "_record"))
        centroid = np.abs(merged["centroid"] - merged["centroid_record"]).max()
        same = (move == 0).all()
        moved = [f"{gene} (row {row + 1}, record {record_row + 1})"
                 for row, (gene, record_row) in enumerate(zip(run["gene"], record_rows)) if row != record_row]
        report("3 figures", species, "MI heatmap rows (centroid order)",
               "identical" if same else "same genes, other order",
               f"{len(run):,} genes; {int((move == 0).sum()):,} in the same row, largest move {int(move.max())} rows, "
               f"Spearman rho of the rows {rho:.6f}; largest centroid difference {centroid:.2g}"
               + (f"; moved: {', '.join(moved)}" if 0 < len(moved) <= 10 else ""))
    for what, run_name, record_name in FIGURES:
        if (figures / run_name).exists():
            result, detail = compare_svg(figures / run_name, package / record_name)
            report("3 figures", "", what, result, detail)

# ---- summary ---------------------------------------------------------------------------------------------------------
table = pd.DataFrame(rows, columns=["step", "species", "what", "result", "detail"])
table.to_csv(args.out, sep="\t", index=False)
print(f"\nRun {args.run} against the heatmaps of record; one row per comparison in {args.out}")
for step in ["1 inputs", "2 monocle", "3 figures"]:
    done = table[table["step"] == step]
    if done.empty:
        print(f"{step}: not run yet")
    else:
        counts = done["result"].value_counts()
        print(f"{step}: " + ", ".join(f"{count} {result}" for result, count in counts.items()))
