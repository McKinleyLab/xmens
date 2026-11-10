"""Step 4: concatenate one species' per-sample QC'd objects into a single object.

Keeps ALL cells (placed + unplaced). Harmonizes obs to use 'sample' (tolerant of
either 'sample' or legacy 'sample_id' on the inputs), adds a clean 'cell_barcode'
with the trailing -1 suffix stripped (for joining to collaborator objects on
sample+barcode), and TRIMS obs to a meaningful whitelist (drops redundant /
constant / vestigial columns - see KEEP_OBS). Preserves layers['counts'] and
obsm['spatial']. Inputs are concatenated in the order given (samples.tsv order in 4_concat.sbatch).
Same code as xmens_final/workflow/scripts/ingest_qc/concat.py (commit 59e00dd), with the Snakemake inputs
replaced by the command-line arguments below.
"""
import argparse
import os
import sys

import anndata as ad
import scanpy as sc

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--species", required=True, choices=["mouse", "human"])
parser.add_argument("--out", required=True, help="output .h5ad; must not exist yet")
parser.add_argument("--log", required=True)
parser.add_argument("inputs", nargs="+", help="per-sample QC'd .h5ad files, in samples.tsv order")
args = parser.parse_args()

inputs = args.inputs
out = args.out
species = args.species
logfile = args.log

if os.path.exists(out):
    sys.exit(f"{out} exists; not overwriting it")
missing = [path for path in inputs if not os.path.exists(path)]
if missing:
    sys.exit(f"missing per-sample objects: {missing}")

# Whitelist of obs columns to keep on the concatenated object. Everything else
# is dropped: n_counts (== total_counts), modality/coords_id (constant/derivable
# from sample), the *_outlier flags (all False post-QC), and log1p_* duplicates.
# 'cell_barcode' is added by this script (below) and always kept.
KEEP_OBS = [
    "sample", "timepoint", "species",
    "x_um_dbscan", "y_um_dbscan",
    "total_counts", "n_genes_by_counts",
    "pct_counts_mt", "pct_counts_ribo",
    "doubletdetection_predicted_doublet", "doubletdetection_doublet_score",
]


def strip_suffix(barcodes):
    """Strip a trailing -<digits> (the 10x GEM-well suffix, e.g. -1)."""
    return barcodes.str.replace(r"-\d+$", "", regex=True)


with open(logfile, "w") as log:
    def say(msg):
        print(f"[concat:{species}] {msg}", file=log, flush=True)

    say(f"concatenating {len(inputs)} samples")

    objs = {}
    for path in inputs:
        a = sc.read_h5ad(path)
        # standardize the sample key: prefer 'sample', fall back to legacy 'sample_id'
        if "sample" not in a.obs.columns:
            if "sample_id" in a.obs.columns:
                a.obs["sample"] = a.obs["sample_id"].astype(str)
            else:
                raise KeyError(f"{path}: no 'sample' or 'sample_id' in obs")
        # drop the legacy column so the concatenated object is clean
        if "sample_id" in a.obs.columns:
            del a.obs["sample_id"]
        sample = str(a.obs["sample"].iloc[0])
        # clean barcode (suffix stripped) for joining to collaborator objects
        a.obs["cell_barcode"] = strip_suffix(a.obs_names.to_series()).values
        say(f"  {sample}: {a.n_obs} cells")
        objs[sample] = a

    # concat keeps all cells; index made unique by appending the sample key so
    # stripped barcodes (unique only within a sample) don't collide.
    adata = ad.concat(
        objs,
        join="outer",
        label="sample_concat",     # records which input each cell came from
        index_unique="-",
        merge="same",              # keep var/obsm columns identical across inputs
    )
    # 'sample' is the authoritative grouping; 'sample_concat' duplicates it - drop.
    if "sample_concat" in adata.obs.columns:
        del adata.obs["sample_concat"]

    adata.obs["species"] = species

    # --- trim obs to the meaningful whitelist (+ cell_barcode, always kept) ---
    keep = [c for c in KEEP_OBS if c in adata.obs.columns]
    keep_set = set(keep) | {"cell_barcode"}
    dropped = [c for c in adata.obs.columns if c not in keep_set]
    adata.obs = adata.obs[[c for c in adata.obs.columns if c in keep_set]].copy()
    say(f"trimmed obs: kept {sorted(keep_set & set(adata.obs.columns))}")
    say(f"trimmed obs: dropped {dropped}")

    say(f"concatenated: {adata.n_obs} cells x {adata.n_vars} genes")
    say(f"obs columns: {list(adata.obs.columns)}")
    say(f"layers: {list(adata.layers.keys())} | obsm: {list(adata.obsm.keys())}")

    # sanity: counts layer + spatial present (scVI + spatial steps need them)
    if "counts" not in adata.layers:
        say("WARNING: no 'counts' layer after concat - scVI needs raw counts")
    if "spatial" not in adata.obsm:
        say("WARNING: no obsm['spatial'] after concat")

    n_placed = int((~adata.obs["x_um_dbscan"].isna()
                    & ~adata.obs["y_um_dbscan"].isna()).sum())
    say(f"placed cells: {n_placed}/{adata.n_obs} "
        f"({100 * n_placed / adata.n_obs:.1f}%)")

    # Written under a temporary name and renamed when complete, so an interrupted run never leaves a file that
    # looks finished.
    partial = os.path.join(os.path.dirname(out) or ".", "." + os.path.basename(out) + ".partial.h5ad")
    adata.write_h5ad(partial)
    os.replace(partial, out)
    say(f"wrote -> {out}")
