"""Compare a run of this pipeline with the SAMap run of record, for each step that has finished.

A step is compared once the run has its outputs, so this can run at any point:
  0 homology    the proteome chunks byte for byte, the phmmer hit tables field by field (their comment lines, which
                hold paths and dates, aside) and the two homology maps, against the run of record (24 Apr 2026);
                "incomplete" while fewer than the 32 chunks or tables of a direction exist
  1 inputs      counts, nuclei and genes, against the inputs of record (inputs_raw/, 1 Sep 2026)
  2 fit         alignment scores and the other numbers of eval.json, the joint graph (edges shared, weights) and the
                joint UMAP, against the fit of record (runs_raw/vanilla_raw)
  3 transfer    Samap_Final.h5ad against the deposited one: obs column by column (for the transferred pseudotime also
                the largest difference and the correlation), every other dataset by checksum
  4 statistics  the per-nucleus tables, graph block, edge shares and permutation control, value by value, against the
                3 Nov 2026 set (record/4_5_figures/results_final/samap_figures)
  5 figures     against the 3 Nov 2026 set: every PNG pixel by pixel, every SVG apart from its drawing date and
                element ids, the CSV value by value
  numbers       the numbers of record/numbers.csv that the run computes (the others are counts drawn in the figures,
                which step 5 compares file by file, and run times)
Inputs the run took from the run of record (start-from variables) are "not compared", as are numbers read unchanged
from such an input.
Writes one row per comparison to --out (CSV) and prints a summary.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem
from PIL import Image
from scipy.stats import pearsonr, spearmanr

NH = "/n/lab_storage/eddy_lab/users/nhilgert"
NS = "/n/netscratch/eddy_lab/Everyone/nhilgert"
FIT = f"{NH}/code/samap-scvi/full_fibdec/runs_raw/vanilla_raw"

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--run", type=Path, required=True, help="the run folder ($OUT)")
parser.add_argument("--phmmer-hits", type=Path, required=True, help="the phmmer hit tables used ($PHMMER_HITS)")
parser.add_argument("--maps", type=Path, required=True, help="the homology maps used ($MAPS, holding hsmm/)")
parser.add_argument("--inputs", type=Path, required=True, help="the SAMap inputs used ($SAMAP_INPUTS)")
parser.add_argument("--joint", type=Path, required=True, help="the joint SAMap object used ($SAMAP_JOINT)")
parser.add_argument("--samap-final", type=Path, required=True, help="the Samap_Final.h5ad used ($SAMAP_FINAL)")
parser.add_argument("--record", type=Path, required=True, help="the record/ folder of this pipeline")
parser.add_argument("--out", type=Path, required=True, help="output CSV; must not exist yet")
parser.add_argument("--record-proteomes", type=Path, default=Path(f"{NH}/data/proteomes"),
                    help="the phmmer run of record: chunks/{hs,mm}/ and hm/{hu_vs_mu,mu_vs_hu}/")
parser.add_argument("--record-maps", type=Path,
                    default=Path(f"{NH}/decidualization/revisions/samap_integrate/homology/hsmm_hmmer"))
parser.add_argument("--record-inputs", type=Path, default=Path(f"{NH}/code/samap-scvi/full_fibdec/inputs_raw"))
parser.add_argument("--record-joint", type=Path, default=Path(f"{FIT}/samap_transferred.h5ad"),
                    help="the joint export of the fit of record")
parser.add_argument("--record-samap-final", type=Path, default=Path(f"{NS}/final_objects/Samap_Final.h5ad"),
                    help="the deposited Samap_Final.h5ad")
args = parser.parse_args()
if args.out.exists():
    raise SystemExit(f"{args.out} exists; not overwriting")

RESULTS = args.record / "4_5_figures" / "results_final"
rows = []


def add(step, item, result, detail="", run="", record=""):
    """One row of the comparison: result is identical, same, differs, incomplete (all present files identical, but
    not all present yet), value (a measurement), missing in run, missing in record or not compared; run and record
    hold the two values where there are two."""
    rows.append({"step": step, "item": item, "result": result, "run": run, "record": record, "detail": detail})


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def svg_normalized(path):
    """The SVG's text without its drawing date, and with the element ids matplotlib draws anew on every run renumbered
    in order of appearance. Only tokens the file defines as ids are renumbered, so colors and paths stay as drawn."""
    text = re.sub(r"<dc:date>[^<]*</dc:date>", "<dc:date/>", Path(path).read_text())
    defined = set(re.findall(r'\bid="([^"]+)"', text))
    names = {}

    def rename(match):
        token = match.group(2)
        if token not in defined:
            return match.group(0)
        return match.group(1) + names.setdefault(token, f"id{len(names)}")

    return re.sub(r'(\bid="|#)([A-Za-z_][\w.-]*)', rename, text)


def same_file(a, b):
    return Path(a).resolve() == Path(b).resolve()


def read_h5ad_parts(path, parts):
    """The named parts of an h5ad file (e.g. "obs", "X", "obsp/connectivities", "obsm/X_umap"), read with h5py."""
    with h5py.File(path, "r") as f:
        return {part: read_elem(f[part]) for part in parts}


def compare_matrices(a, b):
    """'identical', or how two sparse matrices of the same shape differ."""
    a, b = sp.csr_matrix(a, dtype=np.float64), sp.csr_matrix(b, dtype=np.float64)
    if a.shape != b.shape:
        return "differs", f"shapes {a.shape} and {b.shape}"
    difference = a - b
    difference.eliminate_zeros()
    if difference.nnz == 0:
        return "identical", ""
    return "differs", f"{difference.nnz} entries differ, largest by {abs(difference).max():.3g}"


# Columns of the record's tables that this pipeline leaves out on purpose: the embedding_outlier flags, dropped on
# 6 Nov 2026 (no panel has hidden the flagged nuclei since 1 Nov).
DROPPED_COLUMNS = {"embedding_outlier"}


def compare_columns(run, record):
    """'identical', or which columns of two tables (same rows) differ and by how much. Columns the pipeline drops
    on purpose are named in the detail but do not count as a difference."""
    dropped = [c for c in record.columns if c not in run.columns and c in DROPPED_COLUMNS]
    record = record.drop(columns=dropped)
    if list(run.columns) != list(record.columns):
        return "differs", f"columns {list(run.columns)} and {list(record.columns)}"
    if not run.index.equals(record.index):
        return "differs", "rows differ in names or order"
    notes = []
    for column in run.columns:
        a, b = run[column], record[column]
        if pd.api.types.is_float_dtype(a) and pd.api.types.is_float_dtype(b):
            a, b = a.to_numpy(), b.to_numpy()
            if not np.array_equal(np.isnan(a), np.isnan(b)):
                notes.append(f"{column}: missing values differ")
            elif not np.array_equal(a, b, equal_nan=True):
                notes.append(f"{column}: largest difference {np.nanmax(np.abs(a - b)):.3g}")
        elif not a.astype(str).equals(b.astype(str)):
            notes.append(f"{column}: {int((a.astype(str) != b.astype(str)).sum())} values differ")
        categorical_a = isinstance(a.dtype, pd.CategoricalDtype)
        categorical_b = isinstance(b.dtype, pd.CategoricalDtype)
        if categorical_a != categorical_b:
            notes.append(f"{column}: categorical in one table only")
        elif categorical_a and list(a.cat.categories.astype(str)) != list(b.cat.categories.astype(str)):
            notes.append(f"{column}: categories differ (or their order)")
    if notes:
        return "differs", "; ".join(notes)
    return "identical", f"{', '.join(dropped)} left out on purpose" if dropped else ""


def inventory(path, skip=("obs",)):
    """sha256 of every dataset and the attributes of every group outside `skip`, as slim_samap_obs.py checked them."""
    items = {}
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            if any(name == s or name.startswith(s + "/") for s in skip):
                return
            attrs = sorted((k, str(v)) for k, v in obj.attrs.items())
            if isinstance(obj, h5py.Dataset):
                digest = hashlib.sha256()
                if obj.shape == ():
                    value = obj[()]
                    digest.update(str(value).encode() if obj.dtype.kind == "O" else np.asarray(value).tobytes())
                else:
                    for start in range(0, obj.shape[0], 1 << 22):
                        block = obj[start:start + (1 << 22)]
                        digest.update("\x00".join(map(str, block)).encode() if obj.dtype.kind == "O" else block.tobytes())
                items[name] = (str(obj.dtype), obj.shape, attrs, digest.hexdigest())
            else:
                items[name] = ("group", attrs)
        f.visititems(visit)
    return items


def hit_digests(path):
    """sha256 of a phmmer --tblout file's hit lines (its comment lines, which hold paths, options and the date,
    left out): as written, and field by field, ignoring the padding phmmer puts between columns."""
    raw, fields = hashlib.sha256(), hashlib.sha256()
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                raw.update(line)
                fields.update(b" ".join(line.split()) + b"\n")
    return raw.hexdigest(), fields.hexdigest()


def map_lines(path):
    """The lines of a homology map, as a set."""
    with open(path) as f:
        return set(f)


# ---- step 0: homology ------------------------------------------------------------------------------------------
N_CHUNKS = 32


def step0_result(same, compared, found):
    """differs if a file differs; incomplete if all found are the same but fewer than N_CHUNKS were found."""
    if same < compared or compared < found:
        return "differs"
    return "identical" if found == N_CHUNKS else "incomplete"


chunks = args.run / "0_homology" / "chunks"
if chunks.exists():
    for species in ("hs", "mm"):
        found = sorted((chunks / species).glob("*.fa.*"))
        run_chunks = [c for c in found if not c.name.endswith(".partial")]  # .partial: being written, or left by a kill
        record_dir = args.record_proteomes / "chunks" / species
        pairs = [(c, record_dir / c.name) for c in run_chunks if (record_dir / c.name).exists()]
        same = sum(sha256(c) == sha256(r) for c, r in pairs)
        detail = f"{same} of {len(run_chunks)} chunks byte for byte"
        if len(pairs) < len(run_chunks):
            detail += f"; {len(run_chunks) - len(pairs)} without a chunk of record"
        if len(found) > len(run_chunks):
            detail += f"; {len(found) - len(run_chunks)} .partial left out"
        add("0 homology", f"{species} proteome chunks", step0_result(same, len(pairs), len(run_chunks)), detail,
            run=len(run_chunks), record=N_CHUNKS)
if same_file(args.phmmer_hits, args.record_proteomes / "hm"):
    add("0 homology", "phmmer hit tables", "not compared", "the run used the hit tables of record",
        run=str(args.phmmer_hits))
elif args.phmmer_hits.exists():
    for direction in ("hu_vs_mu", "mu_vs_hu"):
        tables = sorted((args.phmmer_hits / direction).glob("*.tbl"))  # not *.tbl.partial
        record_dir = args.record_proteomes / "hm" / direction
        digests = [(hit_digests(t), hit_digests(record_dir / t.name)) for t in tables if (record_dir / t.name).exists()]
        same_fields = sum(run[1] == record[1] for run, record in digests)
        same_bytes = sum(run[0] == record[0] for run, record in digests)
        detail = (f"{same_fields} of {len(tables)} tables with the same hits in every field, "
                  f"{same_bytes} also byte for byte")
        if len(digests) < len(tables):
            detail += f"; {len(tables) - len(digests)} without a table of record"
        add("0 homology", f"phmmer hit tables {direction}", step0_result(same_fields, len(digests), len(tables)),
            detail, run=len(tables), record=N_CHUNKS)
if same_file(args.maps, args.record_maps):
    add("0 homology", "homology maps", "not compared", "the run used the maps of record", run=str(args.maps))
else:
    for name in ("hs_to_mm.txt", "mm_to_hs.txt"):
        run_map, record_map = args.maps / "hsmm" / name, args.record_maps / "hsmm" / name
        if not run_map.exists():
            continue
        if sha256(run_map) == sha256(record_map):
            add("0 homology", f"map {name}", "identical", "byte for byte")
            continue
        a, b = map_lines(run_map), map_lines(record_map)
        detail = f"{len(a - b)} lines only in the run, {len(b - a)} only in the record"
        add("0 homology", f"map {name}", "identical" if a == b else "differs",
            "same lines in another order" if a == b else detail, run=len(a), record=len(b))

# ---- step 1: inputs --------------------------------------------------------------------------------------------
if same_file(args.inputs, args.record_inputs):
    add("1 inputs", "SAMap inputs", "not compared", "the run used the inputs of record", run=str(args.inputs))
else:
    for code in ("hs", "mm"):
        run_file = args.inputs / f"{code}_fib_dec_all_raw.h5ad"
        record_file = args.record_inputs / f"{code}_fib_dec_all_raw.h5ad"
        if not run_file.exists():
            continue
        run_parts = read_h5ad_parts(run_file, ["X", "obs", "var"])
        record_parts = read_h5ad_parts(record_file, ["X", "obs", "var"])
        add("1 inputs", f"{code} counts", *compare_matrices(run_parts["X"], record_parts["X"]))
        same_nuclei = run_parts["obs"].index.astype(str).equals(record_parts["obs"].index.astype(str))
        same_genes = run_parts["var"].index.astype(str).equals(record_parts["var"].index.astype(str))
        add("1 inputs", f"{code} nuclei, in order", "identical" if same_nuclei else "differs",
            run=len(run_parts["obs"]), record=len(record_parts["obs"]))
        add("1 inputs", f"{code} genes, in order", "identical" if same_genes else "differs",
            run=len(run_parts["var"]), record=len(record_parts["var"]))

# ---- step 2: fit ------------------------------------------------------------------------------------------------
# Any object in the folder of the fit of record (the joint export, or the transfer objects made from it) is that fit.
uses_record_fit = args.joint.resolve().parent == args.record_joint.resolve().parent
if uses_record_fit:
    add("2 fit", "joint object", "not compared", "the run used the fit of record", run=str(args.joint))
elif args.joint.exists():
    parts = ["obs", "obsp/connectivities", "obsm/X_umap"]
    run_parts, record_parts = read_h5ad_parts(args.joint, parts), read_h5ad_parts(args.record_joint, parts)
    run_names = run_parts["obs"].index.astype(str)
    record_names = record_parts["obs"].index.astype(str)
    add("2 fit", "nuclei, in order", "identical" if run_names.equals(record_names) else "differs",
        run=len(run_names), record=len(record_names))
    if set(run_names) == set(record_names):
        order = run_names.get_indexer(record_names)
        A = sp.csr_matrix(record_parts["obsp/connectivities"])
        B = sp.csr_matrix(run_parts["obsp/connectivities"])[order][:, order]
        add("2 fit", "joint graph (connectivities)", *compare_matrices(B, A), run=B.nnz, record=A.nnz)
        pattern_a, pattern_b = A.astype(bool).astype(np.int8), B.astype(bool).astype(np.int8)
        shared = pattern_a.multiply(pattern_b).nnz
        add("2 fit", "joint graph: edges shared (Jaccard)", "value", f"{shared} of {A.nnz} record edges",
            run=round(shared / (A.nnz + B.nnz - shared), 4), record=1.0)

        # Human x mouse block: how many of each human nucleus's mouse neighbors are the same.
        species = record_parts["obs"]["species"].astype(str).to_numpy()
        human, mouse = species == "hs", species == "mm"
        block_a = pattern_a[human][:, mouse].tocsr()
        block_b = pattern_b[human][:, mouse].tocsr()
        both = np.asarray(block_a.multiply(block_b).sum(axis=1)).ravel()
        either = np.asarray(block_a.sum(axis=1)).ravel() + np.asarray(block_b.sum(axis=1)).ravel() - both
        jaccard = np.divide(both, either, out=np.ones_like(both, dtype=float), where=either > 0)
        add("2 fit", "mouse neighbors of each human nucleus (mean Jaccard)", "value",
            f"identical for {int((jaccard == 1).sum())} of {len(jaccard)} human nuclei",
            run=round(float(jaccard.mean()), 4), record=1.0)

        umap_a = np.asarray(record_parts["obsm/X_umap"], dtype=float)
        umap_b = np.asarray(run_parts["obsm/X_umap"], dtype=float)[order]
        same_umap = np.array_equal(umap_a, umap_b)
        add("2 fit", "joint UMAP", "identical" if same_umap else "differs",
            "" if same_umap else f"largest coordinate difference {np.abs(umap_a - umap_b).max():.3g}")
    else:
        add("2 fit", "joint graph and UMAP", "not compared", "the run and the record hold different nuclei")

# ---- step 3: Samap_Final.h5ad --------------------------------------------------------------------------------------
record_final = args.record_samap_final
uses_record_final = same_file(args.samap_final, record_final)
if uses_record_final:
    add("3 transfer", "Samap_Final.h5ad", "not compared", "the run used the deposited object",
        run=str(args.samap_final))
elif args.samap_final.exists() and not record_final.exists():
    add("3 transfer", "Samap_Final.h5ad", "missing in record", f"{record_final} not found")
elif args.samap_final.exists():
    run_obs = read_h5ad_parts(args.samap_final, ["obs"])["obs"]
    record_obs = read_h5ad_parts(record_final, ["obs"])["obs"]
    for column in record_obs.columns:
        if column not in run_obs.columns:
            add("3 transfer", f"obs {column}", "missing in run")
            continue
        add("3 transfer", f"obs {column}", *compare_columns(run_obs[[column]], record_obs[[column]]))
    a = run_obs["transferred_mouse_pseudotime"].to_numpy(dtype=float)
    b = record_obs["transferred_mouse_pseudotime"].reindex(run_obs.index).to_numpy(dtype=float)
    both = np.isfinite(a) & np.isfinite(b)
    add("3 transfer", "transferred pseudotime: Pearson r with the record", "value",
        f"{int(both.sum())} human nuclei with both", run=round(float(pearsonr(a[both], b[both]).statistic), 6),
        record=1.0)
    add("3 transfer", "transferred pseudotime: Spearman rho with the record", "value",
        run=round(float(spearmanr(a[both], b[both]).statistic), 6), record=1.0)
    fibdec_a = run_obs["transferred_mouse_fibdec"].astype(str)
    fibdec_b = record_obs["transferred_mouse_fibdec"].reindex(run_obs.index).astype(str)
    human = run_obs["species"].astype(str) == "hs"
    add("3 transfer", "transferred cell type: same as the record", "value", "human nuclei",
        run=int((fibdec_a == fibdec_b)[human].sum()), record=int(human.sum()))

    # uns/provenance is each run's own note (step 3 adds to it), so it is left out of the checksum comparison.
    skip = ("obs", "uns/provenance")
    run_items, record_items = inventory(args.samap_final, skip), inventory(record_final, skip)
    differing = sorted(k for k in set(run_items) | set(record_items) if run_items.get(k) != record_items.get(k))
    add("3 transfer", "every dataset outside obs and uns/provenance", "identical" if not differing else "differs",
        "; ".join(differing[:12]), run=len(run_items), record=len(record_items))
    add("3 transfer", "uns/provenance", "not compared", "each run writes its own note of where its parts came from")

# ---- step 4: statistics -------------------------------------------------------------------------------------------
statistics = args.run / "4_statistics"
record_statistics = RESULTS / "samap_figures"
if statistics.exists():
    for name, reader in [
        ("tables/joint.parquet", lambda p: pd.read_parquet(p)),
        ("tables/mouse_cells.csv", lambda p: pd.read_csv(p, index_col="cell")),
        ("tables/human_cells.csv", lambda p: pd.read_csv(p, index_col="cell")),
        ("edge_share/samap_compartment_edge_share.csv", lambda p: pd.read_csv(p, index_col=[0, 1])),
        ("edge_share/samap_compartment_edge_share_all.csv", lambda p: pd.read_csv(p)),
        ("ptnull/ptnull_statistics_table.csv", lambda p: pd.read_csv(p)),
        ("ptnull/null_replicates.csv", lambda p: pd.read_csv(p)),
        ("ptnull/example_transfers.csv", lambda p: pd.read_csv(p, index_col="cell")),
    ]:
        run_file, record_file = statistics / name, record_statistics / name
        if not record_file.exists():
            add("4 statistics", name, "missing in record")
        elif not run_file.exists():
            add("4 statistics", name, "missing in run")
        else:
            add("4 statistics", name, *compare_columns(reader(run_file), reader(record_file)))
    if (record_statistics / "graph_block" / "C_hs_mm.npz").exists():
        run_block = sp.load_npz(statistics / "graph_block" / "C_hs_mm.npz")
        record_block = sp.load_npz(record_statistics / "graph_block" / "C_hs_mm.npz")
        add("4 statistics", "graph_block/C_hs_mm.npz", *compare_matrices(run_block, record_block))
        for name in ("hs_names.txt", "mm_names.txt"):
            run_text = (statistics / "graph_block" / name).read_text()
            same = run_text == (record_statistics / "graph_block" / name).read_text()
            add("4 statistics", f"graph_block/{name}", "identical" if same else "differs")
    else:
        add("4 statistics", "graph_block/", "missing in record")

    # Fig. 6C and 6E read positions and timepoints from the human plot table.
    if (record_statistics / "tables" / "human_cells.csv").exists():
        used = ["x_um_qc", "y_um_qc", "timepoint"]
        plot_table = pd.read_parquet(statistics / "tables" / "human_plot_table.parquet", columns=used)
        record_cells = pd.read_csv(record_statistics / "tables" / "human_cells.csv", index_col="cell")[used]
        plot_table["timepoint"] = plot_table["timepoint"].astype(str)
        add("4 statistics", "tables/human_plot_table.parquet (columns Fig. 6C/6E read)",
            *compare_columns(plot_table, record_cells))

# ---- step 5: figures ---------------------------------------------------------------------------------------------
FIGURE_FOLDERS = {"samap_main": "fig6", "samap_supp": "sfig_samap"}
figures = args.run / "5_figures"
if figures.exists():
    for record_folder, run_folder in FIGURE_FOLDERS.items():
        if not (RESULTS / record_folder).is_dir():
            add("5 figures", f"{run_folder}/", "missing in record", f"{RESULTS / record_folder} not found")
            continue
        for record_file in sorted((RESULTS / record_folder).rglob("*")):
            if record_file.suffix not in (".png", ".svg", ".csv"):
                continue
            relative = record_file.relative_to(RESULTS / record_folder)
            run_file = figures / run_folder / relative
            item = f"{run_folder}/{relative}"
            if not run_file.exists():
                add("5 figures", item, "missing in run")
            elif record_file.suffix == ".png":
                a, b = np.asarray(Image.open(run_file)), np.asarray(Image.open(record_file))
                if a.shape != b.shape:
                    add("5 figures", item, "differs", "image sizes differ", run=str(a.shape), record=str(b.shape))
                else:
                    changed = int(np.any(a != b, axis=-1).sum()) if a.ndim == 3 else int((a != b).sum())
                    add("5 figures", item, "identical" if changed == 0 else "differs",
                        "" if changed == 0 else f"{changed} of {a.shape[0] * a.shape[1]} pixels differ")
            elif record_file.suffix == ".svg":
                if sha256(run_file) == sha256(record_file):
                    add("5 figures", item, "identical", "byte for byte")
                elif svg_normalized(run_file) == svg_normalized(record_file):
                    add("5 figures", item, "identical", "apart from the drawing date and element ids")
                else:
                    add("5 figures", item, "differs", "beyond the drawing date and element ids")
            else:
                add("5 figures", item, *compare_columns(pd.read_csv(run_file), pd.read_csv(record_file)))

# ---- numbers --------------------------------------------------------------------------------------------------
def flatten(prefix, value, out):
    if isinstance(value, dict):
        for key, item in value.items():
            flatten(f"{prefix}.{key}", item, out)
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        out[prefix] = float(value)


run_numbers = {}
fit_dir = args.joint.parent
fit_meta = json.loads((fit_dir / "meta.json").read_text()) if (fit_dir / "meta.json").exists() else {}
if uses_record_fit:
    add("numbers", "fit numbers (eval.json, meta.json)", "not compared", "the run used the fit of record")
elif "n_iterations" in fit_meta and (fit_dir / "eval.json").exists():  # written by step 2 of this pipeline
    flatten("fit.eval", json.loads((fit_dir / "eval.json").read_text()), run_numbers)
    run_numbers["fit.meta.mean_alignment"] = fit_meta["mean_alignment"]
    run_numbers["fit.meta.n_iter"] = fit_meta["n_iterations"]
    run_numbers["fit.meta.min_expression"] = fit_meta["preprocess"]["min_expression"]
elif args.joint.exists():
    add("numbers", "fit numbers (eval.json, meta.json)", "not compared",
        f"{fit_dir} holds no eval.json and meta.json of this pipeline's step 2")
if uses_record_final:
    add("numbers", "numbers read from Samap_Final.h5ad", "not compared", "the run used the deposited object")
elif args.samap_final.exists():
    obs = read_h5ad_parts(args.samap_final, ["obs"])["obs"]
    human = obs["species"].astype(str) == "hs"
    run_numbers["final_object.human_with_transferred_pseudotime"] = obs.loc[human, "transferred_mouse_pseudotime"].notna().sum()
    run_numbers["final_object.human_with_transferred_fibdec"] = obs.loc[human, "transferred_mouse_fibdec"].notna().sum()
    for label, count in obs.loc[human, "transferred_mouse_fibdec"].value_counts().items():
        run_numbers[f"final_object.transferred_fibdec_counts.{label}"] = count
    run_numbers["final_object.mouse_rows_with_any_transfer"] = (
        ~human & (obs["transferred_mouse_pseudotime"].notna() | obs["transferred_mouse_fibdec"].notna())).sum()
    run_numbers["fig6c.umap_human_with_value"] = run_numbers["final_object.human_with_transferred_pseudotime"]
    run_numbers["fig6c.umap_human_no_value"] = obs.loc[human, "transferred_mouse_pseudotime"].isna().sum()
transfer_summary = args.run / "3_transfer" / "summary.json"
if transfer_summary.exists() and not uses_record_final:
    summary = json.loads(transfer_summary.read_text())
    run_numbers["label_agreement.median_vote_share"] = summary["cell_type"]["median_vote_share"]
if (statistics / "label_agreement.csv").exists():
    agreement = pd.read_csv(statistics / "label_agreement.csv").set_index("native_cell_type")
    run_numbers["label_agreement.agree"] = agreement.loc["all", "agree"]
    run_numbers["label_agreement.percent_of_labelled"] = round(100 * agreement.loc["all", "agree_share_of_labelled"], 2)
    run_numbers["label_agreement.percent_of_all_human"] = round(100 * agreement.loc["all", "agree_share_of_all"], 2)
    for native in agreement.index.drop("all"):
        for key, column in (("agree", "agree"), ("labelled", "labelled"), ("all", "nuclei")):
            run_numbers[f"label_agreement.by_native_cell_type.{native}.{key}"] = agreement.loc[native, column]
if (statistics / "edge_share" / "samap_compartment_edge_share.csv").exists():
    shares = pd.read_csv(statistics / "edge_share" / "samap_compartment_edge_share.csv", index_col=[0, 1])
    for (direction, source), row in shares.iterrows():
        for target, value in row.items():
            run_numbers[f"edge_share.{direction}.{source}.{target}"] = value
if (statistics / "edge_share" / "samap_compartment_edge_weight.csv").exists():
    weights = pd.read_csv(statistics / "edge_share" / "samap_compartment_edge_weight.csv")
    for _, row in weights.iterrows():
        for column in ("mean_total_weight", "max_total_weight", "n_without_cross_edge", "max_cross_edges"):
            run_numbers[f"edge_weight.{row['direction']}.{row['source_compartment']}.{column}"] = row[column]
if (statistics / "alignment_scores.csv").exists():
    scores = pd.read_csv(statistics / "alignment_scores.csv")
    for _, row in scores.iterrows():
        run_numbers[f"alignment_score.{row['human_compartment']}.{row['mouse_compartment']}"] = row["alignment_score"]
if (statistics / "ptnull" / "ptnull_statistics_table.csv").exists():
    table = pd.read_csv(statistics / "ptnull" / "ptnull_statistics_table.csv").set_index("statistic")
    for stat, row in table.iterrows():
        for column, value in row.items():
            run_numbers[f"ptnull.{stat}.{column}"] = value

record_numbers = pd.to_numeric(pd.read_csv(args.record / "numbers.csv").set_index("key")["value"], errors="coerce")
for key, value in sorted(run_numbers.items()):
    if key not in record_numbers.index:
        continue
    recorded = float(record_numbers[key])
    difference = abs(float(value) - recorded)
    result = "same" if difference <= 1e-9 * max(1.0, abs(recorded)) else "differs"
    add("numbers", key, result, "" if result == "same" else f"difference {difference:.3g}", run=float(value),
        record=recorded)

# ---- write and summarize -------------------------------------------------------------------------------------------
table = pd.DataFrame(rows, columns=["step", "item", "result", "run", "record", "detail"])
table.to_csv(args.out, index=False)
print(f"wrote {args.out}: {len(table)} comparisons\n")
print(table.groupby(["step", "result"]).size().unstack(fill_value=0).to_string())
differs = table[table["result"].isin(["differs", "incomplete", "missing in run", "missing in record"])]
if len(differs):
    print("\nnot identical:")
    print(differs.to_string(index=False, max_colwidth=90))
values = table[table["result"] == "value"]
if len(values):
    print("\nmeasured:")
    print(values.to_string(index=False, max_colwidth=90))
