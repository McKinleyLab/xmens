#!/usr/bin/env python
"""Check every value of Supplementary_Table_QC.xlsx against the files it reports on.

Independent of tables/supplementary_table_qc.ipynb: its own parsers for the Cell Ranger and CellBender outputs,
its own implementation of the QC rules on the raw count matrices (not the QC metrics stored in the objects), and
its own version probe. Every comparison is one row of check_mismatches.csv when it fails; check_report.txt
summarizes per sheet. Exit status 1 if anything failed.

    python check_supplementary_table_qc.py Supplementary_Table_QC.xlsx out_folder
"""
import csv
import datetime
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from anndata.io import read_elem

NH = "/n/lab_storage/eddy_lab/users/nhilgert"
NS = "/n/netscratch/eddy_lab/Everyone/nhilgert"
DATA = "/n/lab_storage/eddy_lab/Lab/mckinley/data"
SAMPLE_SHEET = f"{NH}/xmens_scripts/preprocess/samples.tsv"
QC_CONFIG = f"{NH}/xmens_scripts/preprocess/record/config/config.yaml"
CELLRANGER = NS + "/2026_11_06_cellranger_count/{lane}"
CELLBENDER = NH + "/decidualization/xmens_inputs_regen/{sample}/{sample}_cellbender"
QC_RUN = f"{NS}/2026_11_06_qc_regen"
QC_OBJECTS = {"mouse": f"{NS}/2026_11_06_qc_regen_nick_andy_coords/qc_passed_mouse.h5ad",
              "human": f"{NS}/2026_11_06_qc_regen_nick_andy_coords/qc_passed_human.h5ad"}
ALL_CELLS = {"mouse": f"{DATA}/20261110_mouse_allcells_final.h5ad",
             "human": f"{DATA}/20261107_Human_allcells_final.h5ad"}
FIBDEC_OBJECTS = {"mouse": f"{DATA}/20261110_mouse_only_eFib_dec_final.h5ad",
                  "human": f"{DATA}/20261107_Human_fibdec_final_nick_pdt.h5ad"}
FIBDEC_TYPES = {"eFib", "dec"}
ENVS = {"ingest, concat and this notebook (rsc_uv)": f"{NH}/.uv/envs/rsc_uv",
        "QC and DoubletDetection (cuda11.8_env)": f"{NH}/.conda/envs/cuda11.8_env",
        "CellBender (cellbender_022)": f"{NH}/.conda/envs/cellbender_022"}
USER_BASE = "/n/home12/nhilgert/.local"
SPECIES = ["mouse", "human"]
STEPS = ["UMI floor", "mitochondrial fraction", "count-depth MAD", "detected-gene MAD", "doublet"]

workbook_path, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
out_dir.mkdir(parents=True, exist_ok=True)
book = pd.read_excel(workbook_path, sheet_name=None)
mismatches, counts = [], {}


def check(sheet, row, column, found, expected, tol=0.0):
    """One comparison. found is the workbook's value, expected the recomputed one."""
    counts[sheet] = counts.get(sheet, 0) + 1
    ok = False
    if isinstance(expected, (int, float, np.integer, np.floating)) and not isinstance(expected, bool):
        try:
            ok = abs(float(found) - float(expected)) <= tol or (pd.isna(found) and pd.isna(expected))
        except (TypeError, ValueError):
            ok = False
    else:
        ok = str(found) == str(expected) or str(found).lower() == str(expected).lower()
    if not ok:
        mismatches.append({"sheet": sheet, "row": row, "column": column, "workbook": found, "recomputed": expected})


def read_obs(path):
    with h5py.File(path, "r") as f:
        return read_elem(f["obs"])


def library_of(sample):
    return sample.rsplit("_", 1)[1]


# ---------------------------------------------------------------- sources
sheet = pd.read_csv(SAMPLE_SHEET, sep="\t", dtype=str)
sheet["library"] = sheet["lane"].str.replace("_", "")
config = yaml.safe_load(Path(QC_CONFIG).read_text())
MIN_COUNTS, N_MADS, MAX_PCT_MITO = config["qc"]["min_counts"], config["qc"]["n_mads"], config["qc"]["max_pct_mito"]

all_cells = {s: read_obs(p) for s, p in ALL_CELLS.items()}
qc_objects = {s: read_obs(p) for s, p in QC_OBJECTS.items()}
fibdec_index = {s: read_obs(p).index for s, p in FIBDEC_OBJECTS.items()}
object_sample = {}  # library -> (object sample name, timepoint)
for obs in all_cells.values():
    for name, tp in obs[["sample", "timepoint"]].astype(str).drop_duplicates().itertuples(index=False):
        object_sample[library_of(name)] = (name, tp)
sheet["object_sample"] = sheet["library"].map(lambda lib: object_sample[lib][0])
sheet["object_timepoint"] = sheet["library"].map(lambda lib: object_sample[lib][1])


def cellranger_metrics(lane):
    """metrics_summary.csv parsed by the csv module: percentages as floats, counts as ints."""
    with open(Path(CELLRANGER.format(lane=lane)) / "outs" / "metrics_summary.csv", newline="") as handle:
        header, values = list(csv.reader(handle))[:2]
    out = {}
    for name, value in zip(header, values):
        value = value.replace(",", "")
        out[name] = float(value.rstrip("%")) if value.endswith("%") else int(value)
    return out


def cellranger_settings(lane):
    folder = Path(CELLRANGER.format(lane=lane))
    invocation, cmdline = (folder / "_invocation").read_text(), (folder / "_cmdline").read_text()
    versions = json.loads((folder / "_versions").read_text())
    ref = re.search(r"--transcriptome=(\S+)", cmdline).group(1)
    return {"version": versions["pipelines"].replace("cellranger-", ""),
            "chemistry": re.search(r'chemistry\s*=\s*"([^"]+)"', invocation).group(1),
            "include_introns": re.search(r"include_introns\s*=\s*(\w+)", invocation).group(1),
            "reference": Path(ref).name, "cmdline": cmdline}


def qc_metrics(sample, species):
    """Per-barcode total UMIs, detected genes and mitochondrial % from the raw matrix, plus the obs."""
    raw = ad.read_h5ad(f"{QC_RUN}/raw/{sample}.h5ad")
    X = sp.csr_matrix(raw.X)
    prefix = config["mito_prefix"][species]
    mt = np.asarray(raw.var_names.str.startswith(prefix))
    total = np.asarray(X.sum(axis=1)).ravel()
    genes = np.asarray((X > 0).sum(axis=1)).ravel()
    mito = np.asarray(X[:, mt].sum(axis=1)).ravel() / np.where(total > 0, total, np.nan) * 100
    return pd.DataFrame({"umis": total, "genes": genes, "mito": mito}, index=raw.obs_names), raw.obs, int(mt.sum())


def mad_bounds(values):
    """median +/- N_MADS * median absolute deviation (raw MAD, no scaling)."""
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    return median - N_MADS * mad, median + N_MADS * mad


def qc_calls(metrics, kept):
    """The step that removes each barcode, in the pipeline's order; '' for a kept nucleus. Also the cutoffs."""
    step = pd.Series("", index=metrics.index, dtype=object)
    below = metrics["umis"] < MIN_COUNTS
    step[below] = "UMI floor"
    above = metrics.loc[~below]
    limits = {}
    flags = {"mitochondrial fraction": above["mito"] > MAX_PCT_MITO}
    for column, label in (("umis", "count-depth MAD"), ("genes", "detected-gene MAD")):
        values = np.log10(above[column].to_numpy() + 1)
        low, high = mad_bounds(values)
        flags[label] = pd.Series((values < low) | (values > high), index=above.index)
        limits[column] = (10 ** low - 1, 10 ** high - 1)
    remaining = pd.Series(True, index=above.index)
    for label in STEPS[1:4]:
        step[above.index[remaining & flags[label]]] = label
        remaining &= ~flags[label]
    passed = above.index[remaining]
    not_kept = passed.difference(kept)
    step[not_kept] = "doublet"
    extra = kept.difference(passed)
    return step, limits, extra


per_sample = {}
for sample, species, lane in sheet[["sample_id", "species", "lane"]].itertuples(index=False):
    metrics, raw_obs, n_mt = qc_metrics(sample, species)
    kept = read_obs(f"{QC_RUN}/qc/{sample}.h5ad").index
    step, limits, extra = qc_calls(metrics, kept)
    called_csv = pd.read_csv(f"{CELLBENDER.format(sample=sample)}_cell_barcodes.csv", header=None)[0]
    per_sample[sample] = {"metrics": metrics, "raw_obs": raw_obs, "n_mt": n_mt, "kept": kept, "step": step,
                          "limits": limits, "extra_kept": extra, "cellbender_called": set(called_csv.astype(str)),
                          "cr_metrics": cellranger_metrics(lane), "cr": cellranger_settings(lane)}
    print(f"{sample}: {len(metrics)} called, {len(kept)} kept, {int((step != '').sum())} removed, "
          f"{len(extra)} kept nuclei failing a threshold", flush=True)

# ---------------------------------------------------------------- Samples
samples = book["Samples"].set_index("Library")
for sample, species, lane, library, name, tp, sheet_tp in sheet[
        ["sample_id", "species", "lane", "library", "object_sample", "object_timepoint", "timepoint"]].itertuples(
            index=False):
    row, S = samples.loc[library], per_sample[sample]
    cr, m = S["cr"], S["cr_metrics"]
    called = len(S["metrics"])
    check("Samples", library, "Species", row["Species"], species)
    check("Samples", library, "Cell Ranger run", row["Cell Ranger run"], lane)
    check("Samples", library, "Pipeline sample", row["Pipeline sample"], sample)
    check("Samples", library, "Sample (objects)", row["Sample (objects)"], name)
    check("Samples", library, "Timepoint or phase (objects)", row["Timepoint or phase (objects)"], tp)
    check("Samples", library, "Timepoint (sample sheet)", row["Timepoint (sample sheet)"], sheet_tp)
    check("Samples", library, "Reference", row["Reference"], cr["reference"])
    check("Samples", library, "Cell Ranger version", row["Cell Ranger version"], cr["version"])
    check("Samples", library, "Chemistry setting", row["Chemistry setting"], cr["chemistry"])
    check("Samples", library, "Include introns", row["Include introns"], cr["include_introns"])
    for column in row.index:
        if column.startswith("Cell Ranger: "):
            metric = column[len("Cell Ranger: "):].replace(" (%)", "")
            metric = metric[0].upper() + metric[1:]
            metric = {"Q30 bases in barcode": "Q30 Bases in Barcode", "Q30 bases in RNA read": "Q30 Bases in RNA Read",
                      "Q30 bases in UMI": "Q30 Bases in UMI"}.get(metric, metric)
            key = next((k for k in m if k.lower() == metric.lower()), None)
            check("Samples", library, column, row[column], m[key] if key else f"<no metric {metric}>", tol=1e-9)
    check("Samples", library, "CellBender-called nuclei", row["CellBender-called nuclei"], called)
    check("Samples", library, "CellBender-called nuclei (barcodes csv)", called, len(S["cellbender_called"]))
    check("Samples", library, "CellBender-called nuclei (barcodes identical)",
          set(S["metrics"].index.astype(str)) == S["cellbender_called"], True)
    check("Samples", library, "Reads per CellBender-called nucleus", row["Reads per CellBender-called nucleus"],
          round(m["Number of Reads"] / called))
    check("Samples", library, "QC: nuclei above the UMI floor", row["QC: nuclei above the UMI floor"],
          int((S["metrics"]["umis"] >= MIN_COUNTS).sum()))
    cutoff_columns = (("QC: total UMIs lower cutoff", "umis", 0), ("QC: total UMIs upper cutoff", "umis", 1),
                      ("QC: detected genes lower cutoff", "genes", 0), ("QC: detected genes upper cutoff", "genes", 1))
    for column, key, k in cutoff_columns:
        check("Samples", library, column, row[column], S["limits"][key][k], tol=0.051)
    positioned = int((S["raw_obs"]["x_um_dbscan"].notna() & S["raw_obs"]["y_um_dbscan"].notna()).sum())
    check("Samples", library, "Slide-tags: called nuclei with a position",
          row["Slide-tags: called nuclei with a position"], positioned)
    check("Samples", library, "Slide-tags: positioned (% of called)",
          row["Slide-tags: positioned (% of called)"], round(100 * positioned / called, 1), tol=0.051)
    with_row = row["Slide-tags: called nuclei with a coordinates row"]
    check("Samples", library, "Slide-tags: coordinates row >= positioned", with_row >= positioned, True)
    check("Samples", library, "Nuclei after QC", row["Nuclei after QC"], len(S["kept"]))
    in_qc = qc_objects[species].index[qc_objects[species]["sample"].astype(str) == sample]
    check("Samples", library, "Nuclei after QC (QC'd object)", len(S["kept"]), len(in_qc))
    check("Samples", library, "Nuclei after QC (same barcodes)",
          set(S["kept"].astype(str) + f"-{sample}") == set(in_qc.astype(str)), True)
    check("Samples", library, "Nuclei in the all-cells object", row["Nuclei in the all-cells object"],
          int((all_cells[species]["sample"].astype(str) == name).sum()))
    check("Samples", library, "kept nuclei failing a threshold", len(S["extra_kept"]), 0)

# ---------------------------------------------------------------- Waterfall
STEPS_SHEET = "Nuclei_per_step" if "Nuclei_per_step" in book else "Waterfall"
waterfall = book[STEPS_SHEET].set_index("Sample")
recomputed_rows = {}
for sample, species, library, name, tp in sheet[
        ["sample_id", "species", "library", "object_sample", "object_timepoint"]].itertuples(index=False):
    S = per_sample[sample]
    step, raw_obs, kept = S["step"], S["raw_obs"], S["kept"]
    qc_obs = qc_objects[species]
    in_qc = qc_obs.loc[qc_obs["sample"].astype(str) == sample]
    cells = all_cells[species]
    in_cells = cells.loc[cells["sample"].astype(str) == name]
    in_fibdec = in_cells.loc[in_cells.index.isin(fibdec_index[species])]
    positioned = lambda obs: int((obs["x_um_dbscan"].notna() & obs["y_um_dbscan"].notna()).sum())
    rec = {"Species": species, "Timepoint": tp, "Library": library,
           "CellBender-called nuclei": len(step), "Called nuclei with a Slide-tags position": positioned(raw_obs),
           **{f"Removed: {label}": int((step == label).sum()) for label in STEPS},
           "Nuclei after QC": len(kept), "Post-QC nuclei with a Slide-tags position": positioned(in_qc),
           "Removed at annotation": len(in_qc) - len(in_cells), "Nuclei in the all-cells object": len(in_cells),
           "All-cells nuclei with a Slide-tags position": positioned(in_cells),
           "Not in the fibroblast/decidual object": len(in_cells) - len(in_fibdec),
           "Nuclei in the fibroblast/decidual object": len(in_fibdec),
           "Fibroblast/decidual nuclei with a Slide-tags position": positioned(in_fibdec)}
    recomputed_rows[name] = rec
    for column, value in rec.items():
        if column not in waterfall.columns and column.startswith("Not in the fibroblast"):
            column_in_book = "Not in the fibroblast/decidual object (other cell types)"
        else:
            column_in_book = column
        check(STEPS_SHEET, name, column, waterfall.loc[name, column_in_book], value)
    check(STEPS_SHEET, name, "all-cells nuclei are a subset of the QC'd nuclei",
          bool(in_cells.index.isin(in_qc.index).all()), True)
    check(STEPS_SHEET, name, "fib/dec nuclei carry eFib or dec in the all-cells object",
          set(in_fibdec["cell_type"].astype(str)) <= FIBDEC_TYPES, True)
    check(STEPS_SHEET, name, "every eFib/dec nucleus is in the fib/dec object",
          int(in_cells["cell_type"].astype(str).isin(FIBDEC_TYPES).sum()), len(in_fibdec))
numeric = [c for c in waterfall.columns if c not in ("Species", "Timepoint", "Library")]
recomputed_name = {"Not in the fibroblast/decidual object (other cell types)": "Not in the fibroblast/decidual object"}
for species in SPECIES:
    total = waterfall.loc[f"{species} total"]
    mine = [recomputed_rows[n] for n in sheet.loc[sheet["species"] == species, "object_sample"]]
    for column in numeric:
        check(STEPS_SHEET, f"{species} total", column, total[column],
              sum(r[recomputed_name.get(column, column)] for r in mine))
    check(STEPS_SHEET, f"{species} total", "Nuclei after QC = QC'd object", total["Nuclei after QC"],
          len(qc_objects[species]))
    check(STEPS_SHEET, f"{species} total", "all-cells = object", total["Nuclei in the all-cells object"],
          len(all_cells[species]))
    check(STEPS_SHEET, f"{species} total", "fib/dec = object", total["Nuclei in the fibroblast/decidual object"],
          len(fibdec_index[species]))

# ---------------------------------------------------------------- Removed_barcodes
removed = book["Removed_barcodes"]
check("Removed_barcodes", "all", "rows = removed at QC + at annotation", len(removed),
      sum(r[f"Removed: {s}"] for r in recomputed_rows.values() for s in STEPS)
      + sum(r["Removed at annotation"] for r in recomputed_rows.values()))
check("Removed_barcodes", "all", "no duplicate nucleus", int(removed["Nucleus"].duplicated().sum()), 0)
by_sample = removed.groupby("Library")
name_of = dict(zip(sheet["library"], sheet["object_sample"]))
for sample, species, library, name in sheet[["sample_id", "species", "library", "object_sample"]].itertuples(
        index=False):
    S = per_sample[sample]
    rows = by_sample.get_group(library) if library in by_sample.groups else removed.iloc[0:0]
    metrics, step = S["metrics"], S["step"]
    qc_obs, cells = qc_objects[species], all_cells[species]
    expected_qc = set(step.index[step != ""].astype(str) + f"-{sample}")
    in_qc = set(qc_obs.index[qc_obs["sample"].astype(str) == sample].astype(str))
    expected_annot = in_qc - set(cells.index.astype(str))
    listed_qc = set(rows.loc[rows["Stage"] == "QC", "Nucleus"].astype(str))
    listed_annot = set(rows.loc[rows["Stage"] == "annotation", "Nucleus"].astype(str))
    check("Removed_barcodes", name, "QC-stage nuclei = recomputed set", listed_qc == expected_qc, True)
    check("Removed_barcodes", name, "annotation-stage nuclei = QC'd minus all-cells",
          listed_annot == expected_annot, True)
    bad_step = bad_value = bad_reason = 0
    for r in rows.to_dict("records"):
        barcode = r["Nucleus"][: -len(sample) - 1]
        umis, genes, mito = r["Total UMIs"], r["Detected genes"], r["Mitochondrial %"]
        numbers = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", str(r["Reason"]))]
        if r["Stage"] == "QC":
            if step.get(barcode, "") != r["Step"]:
                bad_step += 1
            if r["Step"] == "UMI floor" and not (numbers[0] == umis and numbers[1] == MIN_COUNTS):
                bad_reason += 1
            if r["Step"] == "mitochondrial fraction" and not (abs(numbers[0] - mito) < 0.011
                                                             and numbers[1] == MAX_PCT_MITO):
                bad_reason += 1
            if r["Step"] in ("count-depth MAD", "detected-gene MAD"):
                lo, hi = S["limits"]["umis" if r["Step"] == "count-depth MAD" else "genes"]
                value = umis if r["Step"] == "count-depth MAD" else genes
                if not (numbers[0] == value and abs(numbers[1] - lo) < 0.06 and abs(numbers[2] - hi) < 0.06
                        and (value < lo or value > hi)):
                    bad_reason += 1
        elif r["Step"] != "not in the all-cells object" or r["Nucleus"] not in in_qc:
            bad_step += 1
        m = metrics.loc[barcode]
        same_mito = abs(float(m["mito"]) - mito) < 0.0011 or (np.isnan(m["mito"]) and pd.isna(mito))
        if not (int(m["umis"]) == umis and int(m["genes"]) == genes and same_mito):
            bad_value += 1
        if r["Species"] != species or r["Sample"] != name or r["Cell barcode"] != barcode.split("-")[0]:
            bad_value += 1
    check("Removed_barcodes", name, "rows with a wrong Step", bad_step, 0)
    check("Removed_barcodes", name, "rows whose UMIs/genes/mito %/names differ from the raw matrix", bad_value, 0)
    check("Removed_barcodes", name, "rows whose Reason text disagrees with the values or cutoffs", bad_reason, 0)

# ---------------------------------------------------------------- Post_QC
post = book["Post_QC"].set_index("Sample")
for sample, species, library, name, tp in sheet[
        ["sample_id", "species", "library", "object_sample", "object_timepoint"]].itertuples(index=False):
    cells = all_cells[species]
    in_cells = cells.loc[cells["sample"].astype(str) == name]
    barcodes = in_cells.index.astype(str).str[: -len(sample) - 1]
    m = per_sample[sample]["metrics"].loc[barcodes]
    row = post.loc[name]
    check("Post_QC", name, "Species", row["Species"], species)
    check("Post_QC", name, "Timepoint", row["Timepoint"], tp)
    check("Post_QC", name, "Library", row["Library"], library)
    check("Post_QC", name, "Nuclei", row["Nuclei"], len(in_cells))
    check("Post_QC", name, "Median UMIs", row["Median UMIs"], m["umis"].median(), tol=0.0051)
    check("Post_QC", name, "Mean UMIs", row["Mean UMIs"], m["umis"].mean(), tol=0.0051)
    check("Post_QC", name, "Median detected genes", row["Median detected genes"], m["genes"].median(), tol=0.0051)
    check("Post_QC", name, "Mean detected genes", row["Mean detected genes"], m["genes"].mean(), tol=0.0051)
    check("Post_QC", name, "Median mitochondrial %", row["Median mitochondrial %"], m["mito"].median(), tol=0.0051)
    check("Post_QC", name, "Mean mitochondrial %", row["Mean mitochondrial %"], m["mito"].mean(), tol=0.0051)
    check("Post_QC", name, "Nuclei with a Slide-tags position", row["Nuclei with a Slide-tags position"],
          int((in_cells["x_um_dbscan"].notna() & in_cells["y_um_dbscan"].notna()).sum()))
    stored = qc_objects[species].loc[in_cells.index, ["total_counts", "n_genes_by_counts", "pct_counts_mt"]]
    check("Post_QC", name, "stored QC metrics = raw matrix (max |diff| UMIs)",
          float(np.abs(stored["total_counts"].to_numpy() - m["umis"].to_numpy()).max()), 0.0, tol=1e-6)
    check("Post_QC", name, "stored QC metrics = raw matrix (max |diff| mito %)",
          float(np.abs(stored["pct_counts_mt"].to_numpy() - m["mito"].to_numpy()).max()), 0.0, tol=1e-6)

# ---------------------------------------------------------------- Post_QC_cell_type / subtype
for sheet_name, column in (("Post_QC_cell_type", "cell_type"), ("Post_QC_cell_subtype", "cell_subtype")):
    if sheet_name not in book:
        continue
    table = book[sheet_name]
    label_columns = ["Cell type"] + (["Cell subtype"] if column == "cell_subtype" else [])
    for species in SPECIES:
        cells = all_cells[species]
        keys = ["cell_type"] + (["cell_subtype"] if column == "cell_subtype" else [])
        cross = pd.crosstab([cells[k].astype(str) for k in keys], cells["sample"].astype(str))
        part = table[table["Species"] == species].set_index(label_columns)
        check(sheet_name, species, "number of labels", len(part), len(cross))
        for label, row in part.iterrows():
            key = label if column == "cell_subtype" else label
            try:
                counts_row = cross.loc[key]
            except KeyError:
                check(sheet_name, f"{species}/{label}", "label exists in the object", False, True)
                continue
            for name in sheet.loc[sheet["species"] == species, "object_sample"]:
                check(sheet_name, f"{species}/{label}", name, row[name], int(counts_row.get(name, 0)))
            check(sheet_name, f"{species}/{label}", "Total", row["Total"], int(counts_row.sum()))
        check(sheet_name, species, "Total column sums to the object", int(part["Total"].sum()), len(cells))

# ---------------------------------------------------------------- Parameters
params = book["Parameters"].set_index("Parameter")["Value"]
cmdlines = set()
for lane in sheet["lane"]:
    c = per_sample[sheet.loc[sheet["lane"] == lane, "sample_id"].iloc[0]]["cr"]["cmdline"]
    for flag in ("id", "transcriptome", "fastqs", "sample", "output-dir"):
        c = re.sub(rf"--{flag}=\S+", f"--{flag}=<{flag}>", c)
    cmdlines.add(c.split(" ", 1)[1].strip())
check("Parameters", "Cell Ranger command", "one command for all lanes", len(cmdlines), 1)
check("Parameters", "Cell Ranger command", "Value", params["Cell Ranger command"], next(iter(cmdlines)))
cb_commands = set()
for sample in sheet["sample_id"]:
    lines = Path(f"{CELLBENDER.format(sample=sample)}.log").read_text().splitlines()
    command = lines[next(i for i, l in enumerate(lines) if l.endswith("Command:")) + 1]
    command = re.sub(r"--input \S+", "--input <raw_feature_bc_matrix.h5>", command)
    cb_commands.add(re.sub(r"--output \S+", "--output <sample>/<sample>", command))
check("Parameters", "CellBender command", "one command for all libraries", len(cb_commands), 1)
check("Parameters", "CellBender command", "Value", params["CellBender command"], next(iter(cb_commands)))
check("Parameters", "UMI floor", "Value", params["UMI floor"], MIN_COUNTS)
check("Parameters", "Mitochondrial fraction", "Value", params["Mitochondrial fraction"], MAX_PCT_MITO)
check("Parameters", "MAD bounds", "Value", params["MAD bounds"], N_MADS)
check("Parameters", "p-value threshold", "Value", params["p-value threshold"], config["doublet"]["p_thresh"])
check("Parameters", "Voter threshold", "Value", params["Voter threshold"], config["doublet"]["voter_thresh"])
dd_lines, coord_lines = set(), set()
for sample in sheet["sample_id"]:
    qc_text = Path(f"{QC_RUN}/logs/01_qc_{sample}.log").read_text()
    dd_lines.add(re.search(r"(DoubletDetection [^\n]+)", qc_text).group(1))
    ingest_text = Path(f"{QC_RUN}/logs/00b_ingest_{sample}.log").read_text()
    coord_lines.add(re.search(r"reading coords: \S+ \(([^)]+)\)", ingest_text).group(1))
check("Parameters", "DoubletDetection", "one setting for all samples", len(dd_lines), 1)
check("Parameters", "DoubletDetection", "Value", params["DoubletDetection"], next(iter(dd_lines)))
if "Slide-tags coordinates" in params:
    check("Parameters", "Slide-tags coordinates", "Value", params["Slide-tags coordinates"], next(iter(coord_lines)))
versions_cr = {per_sample[s]["cr"]["version"] for s in sheet["sample_id"]}
check("Parameters", "Cell Ranger count version", "Value", params["Cell Ranger count version"],
      ", ".join(sorted(versions_cr)))
for species in SPECIES:
    refs = {per_sample[s]["cr"]["reference"] for s in sheet.loc[sheet["species"] == species, "sample_id"]}
    check("Parameters", f"Reference ({species})", "Value", params[f"Reference ({species})"], ", ".join(sorted(refs)))

# ---------------------------------------------------------------- Versions
PROBE = ("import importlib, sys\n"
         "for n in sys.argv[1:]:\n"
         "    try:\n        m = importlib.import_module(n)\n    except Exception:\n        continue\n"
         "    v = getattr(m, '__version__', None)\n"
         "    if v is None:\n        import pkg_resources\n"
         "        dist = {'umap': 'umap-learn', 'sklearn': 'scikit-learn'}.get(n, n)\n"
         "        v = pkg_resources.get_distribution(dist).version\n"
         "    print(n, v, m.__file__)\n")
versions = book["Versions"]
env_rows = versions[versions["Environment"].isin(ENVS)]
for label, env in ENVS.items():
    tools = env_rows.loc[env_rows["Environment"] == label, "Tool"].tolist()
    tools = [t for t in tools if t != "python"]
    result = subprocess.run([f"{env}/bin/python", "-c", PROBE, *tools], capture_output=True, text=True, timeout=900,
                            env=dict(os.environ, PYTHONUSERBASE=USER_BASE, PYTHONDONTWRITEBYTECODE="1"))
    probed = {line.split()[0]: line.split()[1:3] for line in result.stdout.splitlines()}
    py = subprocess.run([f"{env}/bin/python", "-c", "import sys; print(sys.version.split()[0])"], capture_output=True,
                        text=True).stdout.strip()
    for r in env_rows[env_rows["Environment"] == label].to_dict("records"):
        tool = r["Tool"]
        if tool == "python":
            check("Versions", f"{label}/python", "Version", r["Version"], py)
            continue
        if tool not in probed:
            check("Versions", f"{label}/{tool}", "importable", False, True)
            continue
        version, file = probed[tool]
        check("Versions", f"{label}/{tool}", "Version", r["Version"], version)
        where = "user site (HOME/.local)" if file.startswith(USER_BASE) else "environment"
        check("Versions", f"{label}/{tool}", "Imported from", r["Imported from"], where)
        check("Versions", f"{label}/{tool}", "Installed", r["Installed"],
              datetime.date.fromtimestamp(os.path.getmtime(file)).isoformat())

# ---------------------------------------------------------------- README
readme = " ".join(book["README"]["Note"].astype(str))
for species in SPECIES:
    for what, n in (("QC'd", len(qc_objects[species])), ("all-cells", len(all_cells[species])),
                    ("fibroblast/decidual", len(fibdec_index[species]))):
        check("README", species, f"{what} count in the objects line", f"({n:,}" in readme, True)
    for path in (QC_OBJECTS[species], ALL_CELLS[species], FIBDEC_OBJECTS[species]):
        check("README", species, f"names {Path(path).name}", Path(path).name in readme, True)

# ---------------------------------------------------------------- report
report = pd.DataFrame(mismatches)
report.to_csv(out_dir / "check_mismatches.csv", index=False)
lines = [f"Check of {workbook_path} on {datetime.datetime.now():%Y-%m-%d %H:%M}", ""]
lines.append(f"{'sheet':24s} {'checks':>7s} {'mismatches':>10s}")
for name in book:
    n = counts.get(name, 0)
    bad = int((report["sheet"] == name).sum()) if len(report) else 0
    lines.append(f"{name:24s} {n:7d} {bad:10d}")
lines += ["", f"total checks {sum(counts.values())}, mismatches {len(report)}"]
if len(report):
    lines += ["", "mismatches (first 200):", report.head(200).to_string(index=False)]
(out_dir / "check_report.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
sys.exit(1 if len(report) else 0)

