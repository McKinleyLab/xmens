"""Step 3a: build one sample's AnnData from its CellBender counts and Slide-tags positions.

Reads the CellBender-corrected count matrix (.h5, no metadata) and the lane's Slide-tags coords CSV, writes
species/timepoint from the command line into obs (the .h5 carries neither), joins the DBSCAN positions by barcode
and keeps ALL cells (coords are nullable obs columns; unplaced cells stay). Writes one .h5ad.
Same code as xmens_final/workflow/scripts/ingest_qc/ingest.py (commit 59e00dd), with the Snakemake inputs
replaced by the command-line arguments below.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--counts", required=True, help="CellBender filtered .h5 (<sample>_cellbender_filtered.h5)")
parser.add_argument("--coords", required=True, help="the lane's Slide-tags coords CSV (<lane>_coords.csv)")
parser.add_argument("--sample", required=True, help="sample_id in samples.tsv")
parser.add_argument("--species", required=True, choices=["mouse", "human"])
parser.add_argument("--timepoint", required=True)
parser.add_argument("--lane", required=True, help="Slide-tags lane, e.g. mck_1; stored as obs coords_id")
parser.add_argument("--modality", default="slidetags")
parser.add_argument("--x-col", default="dbscan_1", help="coords CSV column stored as x_um_dbscan")
parser.add_argument("--y-col", default="dbscan_2", help="coords CSV column stored as y_um_dbscan")
parser.add_argument("--out", required=True, help="output .h5ad; must not exist yet")
parser.add_argument("--log", required=True)
args = parser.parse_args()

counts_path = args.counts
out_path = args.out
sample = args.sample
species = args.species
timepoint = args.timepoint
modality = args.modality
coords_id = args.lane
coords_path = args.coords
coords_x = args.x_col
coords_y = args.y_col
logfile = args.log

if os.path.exists(out_path):
    sys.exit(f"{out_path} exists; not overwriting it")


with open(logfile, "w") as log:
    def say(msg):
        print(f"[ingest:{sample}] {msg}", file=log, flush=True)

    # --- counts (CellBender .h5) ---
    say(f"reading CellBender counts: {counts_path}")
    adata = sc.read_10x_h5(counts_path)
    adata.var_names_make_unique()   # gene symbols have duplicates
    n0 = adata.n_obs
    say(f"loaded {n0} cells x {adata.n_vars} genes")

    # Preserve corrected counts as a layer so downstream (scVI) has raw input
    # even after any future normalization of X.
    adata.layers["counts"] = adata.X.copy()

    # --- metadata from the sample sheet (authoritative) ---
    adata.obs["sample"] = sample
    adata.obs["species"] = species
    adata.obs["timepoint"] = timepoint
    adata.obs["modality"] = modality
    adata.obs["coords_id"] = coords_id
    for c in ["sample", "species", "timepoint", "modality", "coords_id"]:
        adata.obs[c] = adata.obs[c].astype("category")

    # --- coords (Source 1: barcode-keyed fork export) ---
    say(f"reading coords: {coords_path} (x={coords_x}, y={coords_y})")
    co = pd.read_csv(coords_path)

    if "barcode" not in co.columns:
        raise KeyError(f"coords CSV {coords_path} has no 'barcode' column; got {list(co.columns)}")
    for col in (coords_x, coords_y):
        if col not in co.columns:
            raise KeyError(f"coords CSV {coords_path} missing '{col}'; got {list(co.columns)}")

    co = co.set_index("barcode")

    # Direct barcode join (suffixes already match, verified in provenance). If a
    # future sample mismatches, the overlap count below will reveal it.
    common = adata.obs_names.intersection(co.index)
    say(f"barcode overlap: {len(common)}/{n0} cells have a coords row")
    if len(common) == 0:
        raise ValueError(
            f"no barcode overlap between counts and coords for {sample}. "
            f"counts e.g. {list(adata.obs_names[:3])}; coords e.g. {list(co.index[:3])}"
        )

    # Attach DBSCAN coords as nullable obs columns; unplaced cells stay as NaN.
    x = pd.Series(np.nan, index=adata.obs_names, dtype=float)
    y = pd.Series(np.nan, index=adata.obs_names, dtype=float)
    x.loc[common] = pd.to_numeric(co.loc[common, coords_x], errors="coerce")
    y.loc[common] = pd.to_numeric(co.loc[common, coords_y], errors="coerce")
    adata.obs["x_um_dbscan"] = x.to_numpy()
    adata.obs["y_um_dbscan"] = y.to_numpy()

    n_placed = int((~adata.obs["x_um_dbscan"].isna() & ~adata.obs["y_um_dbscan"].isna()).sum())
    say(f"placed cells (non-null dbscan x AND y): {n_placed}/{n0} "
        f"({100 * n_placed / n0:.1f}%)")

    # Also populate obsm['spatial'] for convenience (NaN for unplaced).
    # Spatial rules subset to placed cells before using this.
    adata.obsm["spatial"] = adata.obs[["x_um_dbscan", "y_um_dbscan"]].to_numpy(dtype=float)

    # --- write (ALL cells kept) ---
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    # Written under a temporary name and renamed when complete, so an interrupted run never leaves a file that
    # looks finished.
    partial_path = os.path.join(os.path.dirname(out_path), "." + os.path.basename(out_path) + ".partial.h5ad")
    adata.write_h5ad(partial_path)
    os.replace(partial_path, out_path)
    say(f"wrote {adata.n_obs} cells x {adata.n_vars} genes -> {out_path}")
