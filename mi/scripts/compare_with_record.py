"""Compare a run of this pipeline with the run of record, step by step, for each step the run has finished.

  1 inputs   genes, pseudotime and expression matrix, against the record's genes.npy, pseudotime.npy and expr.npz
  2 chunks   every chunk the run has: genes, MI and p as written, against the record's results/mi_<chunk>.csv
  3 combine  mi_<species>.csv against the record's mi_all.csv: byte for byte, else the rows as text in any order
  4 overlap  the one-to-one pairs, methods_overlap.json and the percentile sweep, against overlap_of_record/
  5 table    every sheet of the workbook, cell by cell, against the approved workbook of 3 Nov 2026. Not compared, as
             they differ by design: README except the inputs' sha256 (its text was rewritten on 12 Nov 2026), and the
             2025 numbers the record carries beside the live values (overlap columns C-D, labels A2 and A27, and
             sheet checks)
  6 fig5f    the selected genes, the pool concordance table and the per-nucleus pool counts, against
             overlap_of_record/crossspecies_heatmap_selection/, and the selection against the 28 genes drawn in Fig. 5F
Writes one row per comparison to --out (TSV) and prints a summary.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

NH = "/n/lab_storage/eddy_lab/users/nhilgert"
RECORD = f"{NH}/xmens_final/results/mi_dpt"
SPECIES = ("mouse", "human")

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--run", type=Path, required=True, help="the run folder ($OUT)")
parser.add_argument("--record-mi", default=f"{RECORD}/pt_harmonized_{{species}}_nall",
                    help="the record's MI folder of each species; {species} is filled in")
parser.add_argument("--record-overlap", type=Path, default=Path(f"{RECORD}/overlap_of_record"))
parser.add_argument("--record-table", type=Path, required=True, help="the approved workbook (record/5_table/)")
parser.add_argument("--panel-genes", type=Path, required=True, help="the 28 genes drawn in Fig. 5F (record/6_fig5f/)")
parser.add_argument("--out", type=Path, required=True, help="output TSV; must not exist yet")
args = parser.parse_args()
if args.out.exists():
    raise SystemExit(f"{args.out} exists; not overwriting")

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


def same_text_rows(run_csv, record_csv, **read_options):
    """'identical' if the files are equal byte for byte, or hold the same table as text, row for row (a gzip file also
    stores its time stamp); 'same rows, other order' if they hold the same rows as text in another order; otherwise
    'DIFFERS'. Also returns the run table as text."""
    if sha256(run_csv) == sha256(record_csv):
        return "identical", pd.read_csv(run_csv, dtype=str, **read_options)
    run = pd.read_csv(run_csv, dtype=str, keep_default_na=False, **read_options)
    record = pd.read_csv(record_csv, dtype=str, keep_default_na=False, **read_options)
    if run.equals(record):
        return "identical", run
    if list(run.columns) == list(record.columns) and \
            sorted(map(tuple, run.values.tolist())) == sorted(map(tuple, record.values.tolist())):
        return "same rows, other order", run
    return "DIFFERS", run


def record_mi(species):
    return Path(args.record_mi.format(species=species))


# ---- 1 inputs ---------------------------------------------------------------------------------------------------------
for species in SPECIES:
    run = args.run / "inputs" / species
    if not run.exists():
        continue
    genes = (run / "genes.txt").read_text().splitlines()
    record_genes = [str(gene) for gene in np.load(record_mi(species) / "genes.npy", allow_pickle=True)]
    report("1 inputs", species, "genes", "identical" if genes == record_genes else "DIFFERS",
           f"{len(genes):,} genes, record {len(record_genes):,}")
    pt, record_pt = np.load(run / "pseudotime.npy"), np.load(record_mi(species) / "pseudotime.npy")
    same = pt.shape == record_pt.shape and np.array_equal(pt, record_pt)
    report("1 inputs", species, "pseudotime", "identical" if same else "DIFFERS",
           f"{len(pt):,} nuclei, record {len(record_pt):,}")
    x, record_x = np.load(run / "expression.npz")["X"], np.load(record_mi(species) / "expr.npz")["X"]
    if x.shape != record_x.shape:
        report("1 inputs", species, "expression", "DIFFERS", f"shape {x.shape}, record {record_x.shape}")
    else:
        n_differ = int((x != record_x).sum())
        detail = f"{x.shape[0]:,} nuclei x {x.shape[1]:,} genes, float32; {n_differ:,} values differ"
        if n_differ:
            detail += f", largest difference {np.abs(x - record_x).max():.2e}"
        report("1 inputs", species, "expression", "identical" if n_differ == 0 else "DIFFERS", detail)
    del x, record_x

# ---- 2 chunks ---------------------------------------------------------------------------------------------------------
for species in SPECIES:
    chunks = sorted((args.run / "chunks" / species).glob("mi_*.csv"))
    if not chunks:
        continue
    results = {path.name: same_text_rows(path, record_mi(species) / "results" / path.name)[0] for path in chunks}
    n_identical = sum(result == "identical" for result in results.values())
    differing = [name for name, result in results.items() if result != "identical"]
    report("2 chunks", species, f"{len(chunks)} chunks", "identical" if not differing else "DIFFERS",
           f"{n_identical} of {len(chunks)} identical (gene, mi, p as written)"
           + (f"; differing: {', '.join(differing[:10])}" if differing else ""))

# ---- 3 combine --------------------------------------------------------------------------------------------------------
for species in SPECIES:
    run_csv = args.run / f"mi_{species}.csv"
    if not run_csv.exists():
        continue
    result, run = same_text_rows(run_csv, record_mi(species) / "mi_all.csv")
    detail = f"{len(run):,} genes"
    if result == "same rows, other order":
        detail += "; same gene, mi, p, q as text, but genes with equal q in another order (numpy sort on this CPU)"
    report("3 combine", species, "mi table", result, detail)

# ---- 4 overlap --------------------------------------------------------------------------------------------------------
overlap = args.run / "overlap"
if overlap.exists():
    same = sha256(overlap / "ortholog_pairs_one2one.csv") == sha256(args.record_overlap / "ortholog_pairs_one2one.csv")
    report("4 overlap", "both", "one-to-one pairs", "identical" if same else "DIFFERS",
           f"{len(pd.read_csv(overlap / 'ortholog_pairs_one2one.csv')):,} pairs")

    run_json = json.loads((overlap / "methods_overlap.json").read_text())
    record_json = json.loads((args.record_overlap / "methods_overlap.json").read_text())
    run_json.pop("mi_tables", None)
    record_json.pop("arms", None)
    flat = lambda d: pd.json_normalize(d).iloc[0].to_dict()                     # noqa: E731 ({"mouse.top": ...})
    run_values, record_values = flat(run_json), flat(record_json)
    differing = [key for key in record_values if run_values.get(key) != record_values[key]]
    detail = (f"N {run_json['N']:,}, K {run_json['K']:,}, n {run_json['n']:,}, k {run_json['k']:,}, "
              f"log10 p {run_json['log10_p']:.6f}")
    if differing:
        detail += "; differing: " + ", ".join(f"{key} {run_values.get(key)} (record {record_values[key]})"
                                              for key in differing)
    report("4 overlap", "both", "methods_overlap.json", "identical" if not differing else "DIFFERS", detail)

    result, _ = same_text_rows(overlap / "overlap_percentile_sweep.csv", args.record_overlap / "overlap_percentile_sweep.csv")
    detail = "60 percentile cutoffs, 40th to 99th"
    if result == "DIFFERS":
        run_sweep = pd.read_csv(overlap / "overlap_percentile_sweep.csv")
        record_sweep = pd.read_csv(args.record_overlap / "overlap_percentile_sweep.csv")
        relative = (run_sweep - record_sweep).abs() / record_sweep.abs().replace(0, 1)
        detail += "; largest relative difference per column: " + ", ".join(
            f"{column} {relative[column].max():.1e}" for column in relative.columns if relative[column].max() > 0)
    report("4 overlap", "both", "percentile sweep", result, detail)


# ---- 5 table ----------------------------------------------------------------------------------------------------------
# The approved workbook carried the 2025 numbers beside the live values (overlap columns C-D, sheet checks) and labels
# that referred to them (overlap A2, A27); since 12 Nov 2026 the workbook has neither, so they are not compared.
OVERLAP_LABELS_REWORDED = {(2, 1), (27, 1)}


def workbook_cells(path):
    """Every non-empty cell of every sheet (formulas as text), without the cells that differ by design."""
    cells = {}
    for ws in openpyxl.load_workbook(path):
        if ws.title == "checks":
            continue
        sheet = {(cell.row, cell.column): cell.value for row in ws.iter_rows() for cell in row if cell.value is not None}
        if ws.title == "README":
            # Only the inputs' sha256 are compared; the paths, build date and descriptions are not.
            sheet = {(row, column): re.findall(r"sha256 ([0-9a-f]{64})", str(value))
                     for (row, column), value in sheet.items()
                     if column == 2 and str(sheet.get((row, 1))).startswith("Input:")}
        if ws.title == "overlap":
            sheet = {key: value for key, value in sheet.items() if key[1] <= 2 and key not in OVERLAP_LABELS_REWORDED}
        cells[ws.title] = sheet
    return cells


workbook = args.run / "table" / "SuppTable_pseudotime_MI.xlsx"
if workbook.exists():
    run_cells, record_cells = workbook_cells(workbook), workbook_cells(args.record_table)
    report("5 table", "both", "sheets", "identical" if list(run_cells) == list(record_cells) else "DIFFERS",
           ", ".join(run_cells))
    for sheet, record_sheet in record_cells.items():
        run_sheet = run_cells.get(sheet, {})
        differing = sorted(key for key in set(record_sheet) | set(run_sheet) if run_sheet.get(key) != record_sheet.get(key))
        detail = f"{len(record_sheet):,} cells, {len(differing):,} differ"
        if differing:
            detail += ": " + "; ".join(f"{openpyxl.utils.get_column_letter(column)}{row} {run_sheet.get((row, column))!r} "
                                       f"(record {record_sheet.get((row, column))!r})" for row, column in differing[:4])
        report("5 table", "both", f"sheet {sheet}", "identical" if not differing else "DIFFERS", detail)

# ---- 6 fig5f ----------------------------------------------------------------------------------------------------------
fig5f = args.run / "fig5f"
record_fig5f = args.record_overlap / "crossspecies_heatmap_selection"
if fig5f.exists():
    result, selected = same_text_rows(fig5f / "conserved_profiles_selected_genes.csv",
                                      record_fig5f / "conserved_profiles_selected_genes.csv")
    report("6 fig5f", "both", "selected genes", result, f"{len(selected)} genes: {', '.join(selected['human'])}")
    result, pool = same_text_rows(fig5f / "conserved_profiles_pool_concordance.csv",
                                  record_fig5f / "conserved_profiles_pool_concordance.csv")
    report("6 fig5f", "both", "pool concordance", result, f"{len(pool)} pool genes")
    for species in SPECIES:
        result, counts = same_text_rows(fig5f / f"percell_{species}_pool.tsv.gz",
                                        record_fig5f / f"percell_{species}_pool.tsv.gz", sep="\t")
        report("6 fig5f", species, "per-nucleus pool counts", result,
               f"{len(counts):,} nuclei x {counts.shape[1] - 3} genes, compared as text (the gzip files store a time stamp)")
    drawn = pd.read_csv(args.panel_genes)["human"].tolist()
    selected_genes = selected["human"].tolist()
    same = selected_genes == drawn
    report("6 fig5f", "both", "genes drawn in Fig. 5F",
           "identical" if same else ("same genes, other order" if sorted(selected_genes) == sorted(drawn) else "DIFFERS"),
           f"{len(set(selected_genes) & set(drawn))} of the {len(drawn)} drawn genes selected")

# ---- summary ----------------------------------------------------------------------------------------------------------
table = pd.DataFrame(rows, columns=["step", "species", "what", "result", "detail"])
table.to_csv(args.out, sep="\t", index=False)
print(f"\nRun {args.run} against the run of record; one row per comparison in {args.out}")
for step in ["1 inputs", "2 chunks", "3 combine", "4 overlap", "5 table", "6 fig5f"]:
    done = table[table["step"] == step]
    if done.empty:
        print(f"{step}: not run yet")
    else:
        counts = done["result"].value_counts()
        print(f"{step}: " + ", ".join(f"{count} {result}" for result, count in counts.items()))
