"""Step 5: name both sets of Slide-tags positions in the QC'd objects: Nick's and Andy's.

Usage: add_coords.py <run folder with qc_passed_{mouse,human}.h5ad> <new destination folder>

Step 3 attaches the positions of Nick's March 2026 Slide-tags run (<nick-dir>/<lane>_coords.csv, DBSCAN x/y) as obs
x_um_dbscan / y_um_dbscan and obsm['spatial']. Before copying anything, this script checks that those fields equal
the CSVs for every nucleus. It then copies each qc_passed object into the new folder, copies those positions to obs
x_um_dbscan_nick / y_um_dbscan_nick, and adds Andy's DBSCAN positions (<andy-dir>/<lane>/Positions: coords.csv rows
indexed into cb_whitelist.txt) as obs x_um_dbscan_andy / y_um_dbscan_andy and obsm['spatial_andy'], joined by cell
barcode. Unplaced nuclei are NaN. Counts and every existing field are copied unchanged, as checked at the end.
Also writes coordinate_summary.tsv (per lane: nuclei placed by each set, by both, and more than 1 um apart).
"""
import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import argparse
import shutil
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("src", type=Path, help="folder with qc_passed_<species>.h5ad (the run folder)")
parser.add_argument("dst", type=Path, help="new folder for the copies; must not exist yet")
parser.add_argument("--nick-dir", type=Path, required=True, help="Nick's Slide-tags coords CSVs (<lane>_coords.csv)")
parser.add_argument("--andy-dir", type=Path, required=True, help="Andy's andy_raw_data folder (<lane>/Positions/)")
parser.add_argument("--species", nargs="+", default=["mouse", "human"], choices=["mouse", "human"])
args = parser.parse_args()

SRC, DST = args.src, args.dst
NICK, ANDY = args.nick_dir, args.andy_dir
SOURCES = {
    "obs x_um_dbscan_nick, y_um_dbscan_nick (same values: obs x_um_dbscan, y_um_dbscan; obsm spatial)":
        "Nick's March 2026 Slide-tags run: " + str(NICK) + "/<lane>_coords.csv, "
        "columns dbscan_1, dbscan_2, attached by ingest.py",
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


def nick_positions(lane):
    """Nick's DBSCAN x/y for one lane (the CSV the pipeline's ingest reads), indexed by the 16-base barcode."""
    pos = pd.read_csv(NICK / f"{lane}_coords.csv", usecols=["barcode", "dbscan_1", "dbscan_2"])
    pos.index = pos["barcode"].astype(str).str[:16].to_numpy()
    assert pos.index.is_unique, lane
    return pos


def prepare(species, summary):
    """Read the source obs (read-only), check the pipeline positions against Nick's CSVs and add both named sets."""
    with h5py.File(SRC / f"qc_passed_{species}.h5ad", "r") as f:
        obs = read_elem(f["obs"])
    x_andy, y_andy = np.full(len(obs), np.nan), np.full(len(obs), np.nan)
    lanes = ("mck_" + obs["sample"].astype(str).str.extract(r"mck_?(\d+[a-z]?)$", expand=False)).to_numpy()
    for lane in sorted(set(lanes)):
        rows = np.flatnonzero(lanes == lane)
        barcodes = obs["cell_barcode"].astype(str).iloc[rows].str[:16]
        pos = andy_positions(lane)
        x_andy[rows] = pos["x_um_dbscan"].reindex(barcodes).to_numpy(float)
        y_andy[rows] = pos["y_um_dbscan"].reindex(barcodes).to_numpy(float)
        x_run, y_run = obs["x_um_dbscan"].to_numpy(float)[rows], obs["y_um_dbscan"].to_numpy(float)[rows]
        nick = nick_positions(lane).reindex(barcodes)
        assert barcodes.isin(nick.index).all(), f"{species} {lane}: nuclei missing from Nick's CSV"
        for run, column, axis in ((x_run, "dbscan_1", "x"), (y_run, "dbscan_2", "y")):
            csv = nick[column].to_numpy(float)
            assert np.allclose(run, csv, rtol=0, atol=1e-6, equal_nan=True), f"{lane}: {axis} != Nick"
        print(f"{species} {lane}: pipeline x/y = Nick's CSV for {len(rows):,} nuclei", flush=True)
        placed_run = ~np.isnan(x_run) & ~np.isnan(y_run)
        placed_andy = ~np.isnan(x_andy[rows]) & ~np.isnan(y_andy[rows])
        both = placed_run & placed_andy
        shift = np.hypot(x_run[both] - x_andy[rows][both], y_run[both] - y_andy[rows][both])
        summary.append({"species": species, "sample": obs["sample"].astype(str).iloc[rows[0]], "lane": lane,
                        "nuclei": len(rows), "placed, Nick coords": int(placed_run.sum()),
                        "placed, Andy coords": int(placed_andy.sum()), "placed by both": int(both.sum()),
                        "both, >1 um apart": int((shift > 1).sum()),
                        "Nick only": int((placed_run & ~placed_andy).sum()),
                        "Andy only": int((placed_andy & ~placed_run).sum())})
    obs["x_um_dbscan_nick"], obs["y_um_dbscan_nick"] = obs["x_um_dbscan"].to_numpy(), obs["y_um_dbscan"].to_numpy()
    obs["x_um_dbscan_andy"], obs["y_um_dbscan_andy"] = x_andy, y_andy
    return obs


def write(species, obs):
    """Copy the source object and replace its obs; add obsm['spatial_andy'] and uns['coordinate_sources']."""
    src, dst = SRC / f"qc_passed_{species}.h5ad", DST / f"qc_passed_{species}.h5ad"
    shutil.copyfile(src, dst)
    with h5py.File(dst, "r+") as f:
        del f["obs"]
        write_elem(f, "obs", obs)
        spatial_andy = obs[["x_um_dbscan_andy", "y_um_dbscan_andy"]].to_numpy(float)
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


rows = []
prepared = {species: prepare(species, rows) for species in args.species}   # all checks before any copy
DST.mkdir()   # fails if it exists, so nothing is overwritten
for species, obs in prepared.items():
    write(species, obs)
summary = pd.DataFrame(rows)
summary.to_csv(DST / "coordinate_summary.tsv", sep="\t", index=False)
print(summary.to_string(index=False))
print(summary.groupby("species")[summary.columns[3:]].sum().to_string())
