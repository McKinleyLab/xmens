"""Compare a run of this pipeline with the November 2026 run of record, sample by sample, for each finished step.

A step is compared once the run has its outputs, so this can run at any point:
  1 Cell Ranger  raw and filtered .h5 and metrics_summary.csv, byte for byte, against the February 2026 run
  2 CellBender   cells called and their counts, against the July 2026 run behind the regen objects
  3 QC           the numbers in each sample's QC log, against the record's logs
  4 concat       nuclei kept and their counts, against the record's qc_passed objects
  5 positions    Nick's and Andy's positions on the shared nuclei, against the regen objects with both sets
Writes one row per sample to --out (TSV) and prints a summary.
"""
import argparse
import hashlib
import re
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

NH = "/n/lab_storage/eddy_lab/users/nhilgert"
SCRATCH = "/n/netscratch/eddy_lab/Everyone/nhilgert"

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--run", type=Path, required=True, help="the run folder ($OUT)")
parser.add_argument("--samples", type=Path, required=True, help="samples.tsv")
parser.add_argument("--out", type=Path, required=True, help="output TSV; must not exist yet")
parser.add_argument("--record-cellranger", default=f"{NH}/decidualization/revisions/velocity/{{lane}}/cellranger/outs")
parser.add_argument("--record-cellbender", default=f"{NH}/decidualization/xmens_inputs_regen/{{sample}}")
parser.add_argument("--record-qc", type=Path, default=Path(f"{SCRATCH}/2026_11_06_qc_regen"),
                    help="the record's QC folder: logs/01_qc_<sample>.log and qc_passed_<species>.h5ad")
parser.add_argument("--record-coords", type=Path, default=Path(f"{SCRATCH}/2026_11_06_qc_regen_nick_andy_coords"),
                    help="the record's qc_passed_<species>.h5ad with both position sets")
args = parser.parse_args()
if args.out.exists():
    raise SystemExit(f"{args.out} exists; not overwriting")

CELLRANGER_FILES = {
    "Cell Ranger raw .h5 identical": "raw_feature_bc_matrix.h5",
    "Cell Ranger filtered .h5 identical": "filtered_feature_bc_matrix.h5",
    "Cell Ranger metrics identical": "metrics_summary.csv",
}
QC_NUMBERS = {   # name -> pattern of the QC log line that holds it
    "loaded": r"loaded (\d+) cells",
    "pass thresholds": r"after threshold QC: (\d+)/",
    "doublets": r"flagged (\d+) doublets",
    "final": r"final: (\d+) cells",
}
POSITIONS = ["x_um_dbscan_nick", "y_um_dbscan_nick", "x_um_dbscan_andy", "y_um_dbscan_andy"]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def identical_rows(a, b):
    """Number of rows (nuclei) with identical counts in a and b, which have the same rows and genes in order."""
    return int(((a != b).getnnz(axis=1) == 0).sum())


def read_10x_h5(path):
    """Cells x genes counts, cell barcodes and gene IDs of a 10x-format .h5 (as CellBender writes it)."""
    with h5py.File(path, "r") as f:
        matrix = f["matrix"]
        n_genes, n_cells = matrix["shape"][:]
        columns = (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:])
        counts = sp.csc_matrix(columns, shape=(n_genes, n_cells))
        return counts.T.tocsr(), pd.Index(matrix["barcodes"][:].astype(str)), matrix["features/id"][:].astype(str)


def compare_cellranger(lane):
    run = args.run / "cellranger" / lane / "outs"
    if not (run / "metrics_summary.csv").exists():
        return {}
    record = Path(args.record_cellranger.format(lane=lane))
    row = {column: sha256(run / name) == sha256(record / name) for column, name in CELLRANGER_FILES.items()}
    row["reads"] = int(str(pd.read_csv(run / "metrics_summary.csv")["Number of Reads"][0]).replace(",", ""))
    return row


def compare_cellbender(sample):
    run_h5 = args.run / "cellbender" / sample / f"{sample}_cellbender_filtered.h5"
    if not run_h5.exists():
        return {}
    record_h5 = Path(args.record_cellbender.format(sample=sample)) / f"{sample}_cellbender_filtered.h5"
    run_counts, run_cells, run_genes = read_10x_h5(run_h5)
    record_counts, record_cells, record_genes = read_10x_h5(record_h5)
    assert np.array_equal(run_genes, record_genes), f"{sample}: CellBender genes differ from the record"
    both = record_cells.intersection(run_cells)
    return {
        "CellBender cells, record": len(record_cells),
        "CellBender cells, run": len(run_cells),
        "CellBender cells, both": len(both),
        "CellBender nuclei with identical counts": identical_rows(record_counts[record_cells.get_indexer(both)],
                                                                   run_counts[run_cells.get_indexer(both)]),
    }


def qc_numbers(log):
    text = log.read_text()
    matches = {name: re.search(pattern, text) for name, pattern in QC_NUMBERS.items()}
    return {name: int(match.group(1)) if match else None for name, match in matches.items()}


def log_lines(log):
    """The log without its "wrote ->" line, whose path differs between runs."""
    return [line for line in log.read_text().splitlines() if "wrote" not in line]


def compare_qc(sample):
    run_log = args.run / "logs" / f"3_qc_{sample}.log"
    if not run_log.exists() or "final:" not in run_log.read_text():
        return {}
    record_log = args.record_qc / "logs" / f"01_qc_{sample}.log"
    record, run = qc_numbers(record_log), qc_numbers(run_log)
    row = {}
    for name in QC_NUMBERS:
        row[f"QC {name}, record"], row[f"QC {name}, run"] = record[name], run[name]
    row["QC log identical"] = log_lines(record_log) == log_lines(run_log)
    return row


def compare_objects(species):
    """Per sample: nuclei in the record's and the run's qc_passed object, in both, and with identical counts."""
    run_path = args.run / f"qc_passed_{species}.h5ad"
    if not run_path.exists():
        return {}
    record, run = ad.read_h5ad(args.record_qc / f"qc_passed_{species}.h5ad"), ad.read_h5ad(run_path)
    assert record.var_names.equals(run.var_names), f"{species}: genes differ from the record"
    both = record.obs_names.intersection(run.obs_names)
    record_rows, run_rows = record.obs_names.get_indexer(both), run.obs_names.get_indexer(both)
    sample_of_both = record.obs["sample"].astype(str).to_numpy()[record_rows]
    rows = {}
    for sample in record.obs["sample"].astype(str).unique():
        mine = np.flatnonzero(sample_of_both == sample)
        rows[sample] = {
            "nuclei kept, record": int((record.obs["sample"].astype(str) == sample).sum()),
            "nuclei kept, run": int((run.obs["sample"].astype(str) == sample).sum()),
            "nuclei kept, both": len(mine),
            "kept nuclei with identical counts": identical_rows(record.layers["counts"][record_rows[mine]],
                                                                run.layers["counts"][run_rows[mine]]),
        }
    return rows


def read_obs(path):
    with h5py.File(path, "r") as f:
        return read_elem(f["obs"])


def compare_positions(species):
    """Per sample: shared nuclei whose four position columns are identical (both unplaced counts as identical)."""
    run_path = args.run / "with_coords" / f"qc_passed_{species}.h5ad"
    if not run_path.exists():
        return {}
    record, run = read_obs(args.record_coords / f"qc_passed_{species}.h5ad"), read_obs(run_path)
    both = record.index.intersection(run.index)
    a, b = record.loc[both, POSITIONS].to_numpy(float), run.loc[both, POSITIONS].to_numpy(float)
    same = ((a == b) | (np.isnan(a) & np.isnan(b))).all(axis=1)
    sample_of_both = record.loc[both, "sample"].astype(str).to_numpy()
    return {sample: {"nuclei with identical positions": int(same[sample_of_both == sample].sum()),
                     "nuclei compared for positions": int((sample_of_both == sample).sum())}
            for sample in np.unique(sample_of_both)}


samples = pd.read_csv(args.samples, sep="\t")
objects = {species: compare_objects(species) for species in ("mouse", "human")}
positions = {species: compare_positions(species) for species in ("mouse", "human")}
rows = []
for sample, lane, species in samples[["sample_id", "lane", "species"]].itertuples(index=False):
    row = {"sample": sample, "lane": lane, "species": species}
    row.update(compare_cellranger(lane))
    row.update(compare_cellbender(sample))
    row.update(compare_qc(sample))
    row.update(objects[species].get(sample, {}))
    row.update(positions[species].get(sample, {}))
    rows.append(row)
    print(f"compared {sample}", flush=True)
table = pd.DataFrame(rows).convert_dtypes()   # whole numbers stay integers despite missing steps
table.to_csv(args.out, sep="\t", index=False)



def overlap(step, record, run, both, identical):
    """One summary line for a step compared nucleus by nucleus; record and run may be None."""
    if both not in table:
        print(f"{step}: not run yet")
        return
    done = table.dropna(subset=[both])
    n_both, n_same = int(done[both].sum()), int(done[identical].sum())
    sizes = f"record {int(done[record].sum()):,}, run {int(done[run].sum()):,}, " if record else ""
    print(f"{step} ({len(done)} samples): {sizes}both {n_both:,}; identical {n_same:,} of {n_both:,} "
          f"({n_same / max(n_both, 1):.2%})")


def identical_in(done, column):
    return f"{int((done[column] == True).sum())}/{len(done)}"   # noqa: E712 (the column also holds NaN)


print(f"\nRun {args.run} against the November 2026 run of record; one row per sample in {args.out}")
if "Cell Ranger raw .h5 identical" in table:
    done = table.dropna(subset=["Cell Ranger raw .h5 identical"])
    parts = [f"{column.replace('Cell Ranger ', '')} {identical_in(done, column)}" for column in CELLRANGER_FILES]
    print(f"1 Cell Ranger ({len(done)} lanes): " + ", ".join(parts))
else:
    print("1 Cell Ranger: not run yet")
overlap("2 CellBender", "CellBender cells, record", "CellBender cells, run", "CellBender cells, both",
        "CellBender nuclei with identical counts")
if "QC log identical" in table:
    done = table.dropna(subset=["QC log identical"])
    print(f"3 QC ({len(done)} samples): log identical in {identical_in(done, 'QC log identical')}; numbers:")
    for _, row in done.iterrows():
        changes = [f"{name} {row[f'QC {name}, record']:.0f} -> {row[f'QC {name}, run']:.0f}"
                   for name in QC_NUMBERS if row[f"QC {name}, record"] != row[f"QC {name}, run"]]
        print(f"    {row['sample']}: " + ("; ".join(changes) if changes else "same"))
else:
    print("3 QC: not run yet")
overlap("4 concat", "nuclei kept, record", "nuclei kept, run", "nuclei kept, both", "kept nuclei with identical counts")
overlap("5 positions", None, None, "nuclei compared for positions", "nuclei with identical positions")
