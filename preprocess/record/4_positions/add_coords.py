"""Name both sets of Slide-tags positions in the February re-run's QC'd objects: Nick's and Andy's.

The pipeline's ingest attaches the positions of the March 2026 Slide-tags re-run (Slide-tags/coords/<lane>_coords.csv,
DBSCAN x/y) as obs x_um_dbscan / y_um_dbscan and obsm['spatial']. This script copies each qc_passed object into a new
folder, copies the pipeline's positions to obs x_um_dbscan_nick / y_um_dbscan_nick, and adds Andy's DBSCAN positions
(andy_raw_data/<lane>/Positions: coords.csv rows indexed into cb_whitelist.txt) as obs x_um_dbscan_andy /
y_um_dbscan_andy and obsm['spatial_andy'], joined by cell barcode. Unplaced nuclei are NaN. Counts and every existing
field are copied unchanged (x_um_dbscan, y_um_dbscan and obsm['spatial'] stay for the pipeline), as checked at the end.
"""
import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import shutil
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

SRC = Path("/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_qc_regen")
DST = Path("/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_qc_regen_nick_andy_coords")
ANDY = Path("/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data/andy_raw_data")
SOURCES = {
    "obs x_um_dbscan_nick, y_um_dbscan_nick (same values: obs x_um_dbscan, y_um_dbscan; obsm spatial)":
        "March 2026 Slide-tags re-run: /n/lab_storage/eddy_lab/users/nhilgert/Slide-tags/coords/<lane>_coords.csv, "
        "columns dbscan_1, dbscan_2, attached by the pipeline's ingest.py",
    "obs x_um_dbscan_andy, y_um_dbscan_andy (same values: obsm spatial_andy)":
        "Andy's original positioning: " + str(ANDY) + "/<lane>/Positions/coords.csv, columns x_um_dbscan, "
        "y_um_dbscan, rows mapped to barcodes through cb_whitelist.txt (cb_index is 1-based)",
}


def andy_positions(lane):
    """Andy's DBSCAN x/y for one lane, indexed by the 16-base barcode."""
    whitelist = pd.read_csv(ANDY / lane / "Positions" / "cb_whitelist.txt", header=None)[0].astype(str).str[:16]
    pos = pd.read_csv(ANDY / lane / "Positions" / "coords.csv", usecols=["cb_index", "x_um_dbscan", "y_um_dbscan"])
    pos.index = whitelist.iloc[pos["cb_index"].to_numpy() - 1].to_numpy()
    assert pos.index.is_unique, lane
    return pos


def add_coords(species, summary):
    src, dst = SRC / f"qc_passed_{species}.h5ad", DST / f"qc_passed_{species}.h5ad"
    shutil.copyfile(src, dst)
    with h5py.File(dst, "r+") as f:
        obs = read_elem(f["obs"])
        x_andy, y_andy = np.full(len(obs), np.nan), np.full(len(obs), np.nan)
        lanes = ("mck_" + obs["sample"].astype(str).str.extract(r"mck_?(\d+[a-z]?)$", expand=False)).to_numpy()
        for lane in sorted(set(lanes)):
            rows = np.flatnonzero(lanes == lane)
            pos = andy_positions(lane)
            barcodes = obs["cell_barcode"].astype(str).iloc[rows].str[:16]
            x_andy[rows] = pos["x_um_dbscan"].reindex(barcodes).to_numpy(float)
            y_andy[rows] = pos["y_um_dbscan"].reindex(barcodes).to_numpy(float)
            x_run, y_run = obs["x_um_dbscan"].to_numpy(float)[rows], obs["y_um_dbscan"].to_numpy(float)[rows]
            placed_run = ~np.isnan(x_run) & ~np.isnan(y_run)
            placed_andy = ~np.isnan(x_andy[rows]) & ~np.isnan(y_andy[rows])
            both = placed_run & placed_andy
            shift = np.hypot(x_run[both] - x_andy[rows][both], y_run[both] - y_andy[rows][both])
            summary.append({"species": species, "sample": obs["sample"].astype(str).iloc[rows[0]], "lane": lane,
                            "nuclei": len(rows), "placed, re-run coords": int(placed_run.sum()),
                            "placed, Andy coords": int(placed_andy.sum()), "placed by both": int(both.sum()),
                            "both, >1 um apart": int((shift > 1).sum()),
                            "re-run only": int((placed_run & ~placed_andy).sum()),
                            "Andy only": int((placed_andy & ~placed_run).sum())})
        obs["x_um_dbscan_nick"], obs["y_um_dbscan_nick"] = obs["x_um_dbscan"].to_numpy(), obs["y_um_dbscan"].to_numpy()
        obs["x_um_dbscan_andy"], obs["y_um_dbscan_andy"] = x_andy, y_andy
        del f["obs"]
        write_elem(f, "obs", obs)
        spatial_andy = np.column_stack([x_andy, y_andy])
        if "obsm" in f:
            write_elem(f["obsm"], "spatial_andy", spatial_andy)
        else:
            write_elem(f, "obsm", {"spatial_andy": spatial_andy})
        if "uns" in f:
            write_elem(f["uns"], "coordinate_sources", SOURCES)
        else:
            write_elem(f, "uns", {"coordinate_sources": SOURCES})
    check(src, dst)


def check(src, dst):
    """Everything copied unchanged except the new coordinate fields."""
    with h5py.File(src, "r") as a, h5py.File(dst, "r") as b:
        for key in ("X/data", "X/indices", "X/indptr", "layers/counts/data", "layers/counts/indices",
                    "layers/counts/indptr"):
            if key in a:
                assert np.array_equal(a[key][:], b[key][:]), f"{dst.name}: {key} differs"
        obs_a, obs_b = read_elem(a["obs"]), read_elem(b["obs"])
        new_columns = ["x_um_dbscan_nick", "y_um_dbscan_nick", "x_um_dbscan_andy", "y_um_dbscan_andy"]
        pd.testing.assert_frame_equal(obs_a, obs_b.drop(columns=new_columns))
        for axis in ("x", "y"):
            assert np.array_equal(obs_b[f"{axis}_um_dbscan_nick"].to_numpy(), obs_a[f"{axis}_um_dbscan"].to_numpy(),
                                  equal_nan=True), f"{dst.name}: {axis}_um_dbscan_nick"
        obsm_a = set(a["obsm"].keys()) if "obsm" in a else set()
        assert set(b["obsm"].keys()) == obsm_a | {"spatial_andy"}, f"{dst.name}: obsm keys"
    back = ad.read_h5ad(dst, backed="r")
    print(f"{dst.name}: {back.n_obs:,} nuclei x {back.n_vars:,} genes; obsm {list(back.obsm.keys())}; "
          f"uns {list(back.uns.keys())}; checks passed", flush=True)
    back.file.close()


DST.mkdir()   # fails if it exists, so nothing is overwritten
rows = []
for species in ("mouse", "human"):
    add_coords(species, rows)
summary = pd.DataFrame(rows)
summary.to_csv(DST / "coordinate_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(summary.to_string(index=False))
print(summary.groupby("species")[summary.columns[3:]].sum().to_string())
