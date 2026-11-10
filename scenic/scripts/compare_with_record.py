#!/usr/bin/env python
"""Compares a run of scenic with the SCENIC run of record, step by step, for every step that has finished.

Run folder (config.sh OUT): scenic_seeds/<species>/ (steps 1-4), correlations/ (step 5), figures/ (step 6).
Record: --record, the September-October 2026 lane (scenic_seeds/<species>/ with the same layout), and --record-outputs
(record/outputs/ of the pipeline folder: the Fig. 5E table and per-species rho of record, the Fig. 5E and Fig. S28 SVGs
of the manuscript, and the per-run regulon sizes behind Fig. S28; record/SOURCES.tsv lists where each came from).

One row per comparison: step, species, seed, item, result, detail. result is "identical" (same values; "same bytes"
in the detail where the files are byte-identical too), "differs", "not run" (the run has no such output yet) or
"no record". A sketch's counts are not compared separately: they are the pooled loom's counts for the sketch's nuclei,
so identical nuclei and an identical pooled loom make them identical.

Usage (in the pySCENIC environment: the regulons are pickled pySCENIC objects):
  python compare_with_record.py --run <OUT> --record <lane> --record-outputs record/outputs --runs runs.tsv --out <tsv>
"""
import argparse
import collections
import hashlib
import json
import os
import pickle
import re

import h5py
import numpy as np
import pandas as pd

TAG = "sketch15k"
CONS6 = f"consensus_{TAG}_weighted_soft_occ6"
CONS2 = f"consensus_{TAG}_weighted_soft"
AUCELL = f"{CONS6}_pyscenic_cli"
SPECIES = ("mouse", "human")
ROWS = []


def add(step, species, seed, item, result, detail="", quiet=False):
    """One row of the output table; printed to the log unless quiet (the runs not yet made, 2 x 200 rows)."""
    ROWS.append({"step": step, "species": species, "seed": seed, "item": item, "result": result, "detail": detail})
    if not quiet:
        tail = f" ({detail})" if detail else ""
        print(f"step {step} {species} {seed} {item}: {result}{tail}".replace("  ", " "), flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def strs(ds):
    try:
        return np.asarray(ds.asstr()[:], dtype=object)
    except (TypeError, AttributeError):
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in ds[:]], dtype=object)


def h5ad_index(f, group):
    g = f[group]
    return strs(g[g.attrs["_index"]])


def sparse_equal(a, b):
    """Two anndata sparse matrices (h5py groups) or dense datasets hold the same values."""
    import scipy.sparse as sp

    def load(node):
        if isinstance(node, h5py.Dataset):
            return sp.csr_matrix(node[:])
        enc = str(node.attrs.get("encoding-type", "csr_matrix"))
        shape = tuple(int(x) for x in node.attrs["shape"])
        cls = sp.csr_matrix if "csr" in enc else sp.csc_matrix
        return cls((node["data"][:], node["indices"][:], node["indptr"][:]), shape=shape).tocsr()

    ma, mb = load(a), load(b)
    if ma.shape != mb.shape:
        return False, f"shapes {ma.shape} and {mb.shape}"
    n = (ma != mb).nnz
    return n == 0, f"{ma.shape[0]:,} x {ma.shape[1]:,}, {ma.nnz:,} non-zero" if n == 0 else f"{n:,} values differ"


def files_equal(run_file, record_file):
    """(identical?, detail) for two files: same bytes, else None to let the caller compare values."""
    if sha256(run_file) == sha256(record_file):
        return True, "same bytes"
    return None, ""


def frames_equal(a, b):
    """Two data frames with the same labels hold the same values; returns (identical, detail)."""
    if list(a.columns) != list(b.columns) or len(a) != len(b) or not a.index.equals(b.index):
        return False, f"labels differ ({a.shape} vs {b.shape})"
    num = a.select_dtypes("number").columns
    other = [c for c in a.columns if c not in num]
    for c in other:
        if not a[c].astype(str).equals(b[c].astype(str)):
            return False, f"column {c} differs"
    if len(num):
        x, y = a[num].to_numpy(float), b[num].to_numpy(float)
        same = (x == y) | (np.isnan(x) & np.isnan(y))
        if not same.all():
            return False, f"{int((~same).sum()):,} values differ, largest difference {np.nanmax(np.abs(x - y)):.3g}"
    return True, f"{a.shape[0]:,} rows x {a.shape[1]} columns, same values"


def load_regulons(path):
    with open(path, "rb") as fh:
        regs = pickle.load(fh)
    return regs


def regulon_key(r, occurrence=False):
    items = tuple(sorted((g, float(w)) for g, w in r.gene2weight.items()))
    occ = tuple(sorted((g, int(n)) for g, n in r.gene2occurrence.items())) if occurrence else ()
    return (r.name, items, occ)


def compare_regulons(run_file, record_file, occurrence=False):
    a = sorted(regulon_key(r, occurrence) for r in load_regulons(run_file))
    b = sorted(regulon_key(r, occurrence) for r in load_regulons(record_file))
    if a == b:
        return True, f"{len(a)} regulons, same targets and weights" + (" and run counts" if occurrence else "")
    na, nb = {k[0] for k in a}, {k[0] for k in b}
    same_names = len(na & nb)
    differing = sum(1 for x, y in zip(a, b) if x != y)
    return False, f"{len(a)} vs {len(b)} regulons, {same_names} names shared, {differing} differ"


# ---------------------------------------------------------------- step 1: pooled objects
def step1(run, rec, sp):
    rd, cd = f"{run}/scenic_seeds/{sp}", f"{rec}/scenic_seeds/{sp}"
    man = f"{rd}/{sp}_dsc_pooled.manifest.json"
    if not os.path.exists(man):
        add(1, sp, "", "pooled objects", "not run")
        return
    a, b = json.load(open(man)), json.load(open(f"{cd}/{sp}_dsc_pooled.manifest.json"))
    keys = [k for k in ("n_cells_full", "n_genes_full", "gene_filter_min_cells", "n_genes_grn", "n_tfs_in_list",
                        "n_tfs_retained", "raw_source") if k in b]
    diff = {k: (a.get(k), b[k]) for k in keys if a.get(k) != b[k]}
    add(1, sp, "", "manifest", "identical" if not diff else "differs",
        ", ".join(f"{k} {b[k]}" for k in keys) if not diff else f"run vs record: {diff}")

    with h5py.File(f"{rd}/{sp}_dsc_pooled.loom", "r") as A, h5py.File(f"{cd}/{sp}_dsc_pooled.loom", "r") as B:
        ga, gb = strs(A["row_attrs/Gene"]), strs(B["row_attrs/Gene"])
        ca, cb = strs(A["col_attrs/CellID"]), strs(B["col_attrs/CellID"])
        same_g = len(ga) == len(gb) and bool((ga == gb).all())
        same_c = len(ca) == len(cb) and bool((ca == cb).all())
        add(1, sp, "", "pooled loom: genes", "identical" if same_g else "differs", f"{len(ga):,} vs {len(gb):,}")
        add(1, sp, "", "pooled loom: nuclei", "identical" if same_c else "differs", f"{len(ca):,} vs {len(cb):,}")
        if same_g and same_c:
            MA, MB = A["matrix"], B["matrix"]
            n_diff = 0
            for j in range(0, MA.shape[1], 2000):
                n_diff += int((MA[:, j:j + 2000] != MB[:, j:j + 2000]).sum())
            add(1, sp, "", "pooled loom: counts", "identical" if n_diff == 0 else "differs",
                f"{MA.shape[0]:,} genes x {MA.shape[1]:,} nuclei" if n_diff == 0 else f"{n_diff:,} values differ")

    with h5py.File(f"{rd}/{sp}_dsc_pooled.h5ad", "r") as A, h5py.File(f"{cd}/{sp}_dsc_pooled.h5ad", "r") as B:
        same = all(np.array_equal(h5ad_index(A, g), h5ad_index(B, g)) for g in ("obs", "var"))
        add(1, sp, "", "pooled h5ad: nuclei and genes", "identical" if same else "differs")
        if same:
            ok, detail = sparse_equal(A["layers/raw"], B["layers/raw"])
            add(1, sp, "", "pooled h5ad: raw counts", "identical" if ok else "differs", detail)
            sa, sb = A["obsm/X_scVI"][:], B["obsm/X_scVI"][:]
            ok = sa.shape == sb.shape and np.array_equal(sa, sb)
            add(1, sp, "", "pooled h5ad: X_scVI (the sketches' embedding)", "identical" if ok else "differs",
                f"{sa.shape[0]:,} x {sa.shape[1]}")


# ---------------------------------------------------------------- steps 2 and 3: the runs
def step2(run, rec, sp, seed):
    rel_sk = f"scenic_seeds/{sp}/looms_{TAG}/{sp}_dsc_{TAG}_seed{seed}.loom"
    rel_adj = f"scenic_seeds/{sp}/adjacencies_{TAG}/adjacencies_seed{seed}.tsv"
    found = False
    if os.path.exists(f"{run}/{rel_sk}"):
        found = True
        with h5py.File(f"{run}/{rel_sk}", "r") as A, h5py.File(f"{rec}/{rel_sk}", "r") as B:
            ca, cb = strs(A["col_attrs/CellID"]), strs(B["col_attrs/CellID"])
            ga, gb = strs(A["row_attrs/Gene"]), strs(B["row_attrs/Gene"])
        ok = len(ca) == len(cb) and bool((ca == cb).all()) and len(ga) == len(gb) and bool((ga == gb).all())
        add(2, sp, seed, "sketch: nuclei and genes", "identical" if ok else "differs",
            f"{len(ca):,} nuclei, {len(ga):,} genes" if ok else f"{len(set(ca) & set(cb)):,} of {len(cb):,} nuclei shared")
    if os.path.exists(f"{run}/{rel_adj}"):
        found = True
        a = pd.read_csv(f"{run}/{rel_adj}", sep="\t")
        b = pd.read_csv(f"{rec}/{rel_adj}", sep="\t")
        m = a.merge(b, on=["TF", "target"], how="outer", suffixes=("_run", "_record"), indicator=True)
        only_run, only_rec = int((m._merge == "left_only").sum()), int((m._merge == "right_only").sum())
        both = m[m._merge == "both"]
        n_imp = int((both.importance_run != both.importance_record).sum())
        if only_run == only_rec == n_imp == 0:
            add(2, sp, seed, "adjacencies", "identical", f"{len(b):,} links, same importances")
        else:
            biggest = float((both.importance_run - both.importance_record).abs().max()) if len(both) else float("nan")
            add(2, sp, seed, "adjacencies", "differs",
                f"{len(b):,} links in the record; {only_run:,} only in the run, {only_rec:,} only in the record, "
                f"{n_imp:,} shared links with another importance (largest difference {biggest:.3g})")
    if not found:
        add(2, sp, seed, "sketch and adjacencies", "not run", quiet=True)


def step3(run, rec, sp, seed):
    rel = f"scenic_seeds/{sp}/regulons_{TAG}/regulons_seed{seed}.pkl"
    if not os.path.exists(f"{run}/{rel}"):
        add(3, sp, seed, "regulons", "not run", quiet=True)
        return
    ok, detail = compare_regulons(f"{run}/{rel}", f"{rec}/{rel}")
    add(3, sp, seed, "regulons", "identical" if ok else "differs", detail)


# ---------------------------------------------------------------- step 4: consensus and AUCell
def step4(run, rec, sp):
    for folder, floor in ((CONS6, 6), (CONS2, 2)):
        rd, cd = f"{run}/scenic_seeds/{sp}/{folder}", f"{rec}/scenic_seeds/{sp}/{folder}"
        if not os.path.isdir(rd):
            add(4, sp, "", f"consensus, floor {floor}", "not run")
            continue
        ok, detail = files_equal(f"{rd}/recurrence_report.tsv", f"{cd}/recurrence_report.tsv")
        if ok is None:
            ok, detail = frames_equal(pd.read_csv(f"{rd}/recurrence_report.tsv", sep="\t"),
                                      pd.read_csv(f"{cd}/recurrence_report.tsv", sep="\t"))
        else:
            detail += f", {len(pd.read_csv(f'{cd}/recurrence_report.tsv', sep=chr(9))):,} TFs"
        add(4, sp, "", f"floor {floor}: recurrence report", "identical" if ok else "differs", detail)
        ok, detail = compare_regulons(f"{rd}/consensus_regulons.pkl", f"{cd}/consensus_regulons.pkl", occurrence=True)
        add(4, sp, "", f"floor {floor}: consensus regulons", "identical" if ok else "differs", detail)
        for name, sep in (("auc_mtx.csv", ","), ("auc_by_timepoint.tsv", "\t")):
            ok, detail = files_equal(f"{rd}/{name}", f"{cd}/{name}")
            if ok is None:
                ok, detail = frames_equal(pd.read_csv(f"{rd}/{name}", sep=sep, index_col=0),
                                          pd.read_csv(f"{cd}/{name}", sep=sep, index_col=0))
            add(4, sp, "", f"floor {floor}: {name} (comparator AUCell, all genes)", "identical" if ok else "differs",
                detail)
    rd, cd = f"{run}/scenic_seeds/{sp}/{AUCELL}", f"{rec}/scenic_seeds/{sp}/{AUCELL}"
    if not os.path.isdir(rd):
        add(4, sp, "", "AUCell of record", "not run")
        return
    ok, detail = files_equal(f"{rd}/consensus_regulons.yaml", f"{cd}/consensus_regulons.yaml")
    add(4, sp, "", "AUCell: consensus_regulons.yaml", "identical" if ok else "differs", detail or "bytes differ")
    for name in ("auc_mtx_weighted.csv", "auc_mtx_unweighted.csv"):
        ok, detail = files_equal(f"{rd}/{name}", f"{cd}/{name}")
        if ok is None:
            ok, detail = frames_equal(pd.read_csv(f"{rd}/{name}", index_col=0), pd.read_csv(f"{cd}/{name}", index_col=0))
        else:
            with open(f"{cd}/{name}") as fh:
                n_reg = len(fh.readline().rstrip("\n").split(",")) - 1
                n_nuc = sum(1 for _ in fh)
            detail += f", {n_nuc:,} nuclei x {n_reg} regulons"
        label = "the activity of record" if "unweighted" not in name else "for comparison"
        add(4, sp, "", f"AUCell: {name} ({label})", "identical" if ok else "differs", detail)


# ---------------------------------------------------------------- step 5: correlations
def quadrant_counts(t, threshold=0.2):
    """Pairs of the Fig. 5E table with |rho| >= threshold in both species, by direction."""
    h, m = t.rho_human.to_numpy(), t.rho_mouse.to_numpy()
    both = (np.abs(h) >= threshold) & (np.abs(m) >= threshold)
    return {"up in both": int((both & (h > 0) & (m > 0)).sum()), "down in both": int((both & (h < 0) & (m < 0)).sum()),
            "up human, down mouse": int((both & (h > 0) & (m < 0)).sum()),
            "down human, up mouse": int((both & (h < 0) & (m > 0)).sum())}


def step5(run, recout):
    d = f"{run}/correlations"
    if not os.path.isdir(d):
        add(5, "both", "", "correlations", "not run")
        return
    for sp in SPECIES:
        a = pd.read_csv(f"{d}/regulon_rho_pseudotime_{sp}.tsv", sep="\t").set_index("regulon")
        b = pd.read_csv(f"{recout}/regulon_rho_pseudotime_{sp}_fibdec_latent.tsv", sep="\t").set_index("regulon")
        ok, detail = frames_equal(a.sort_index(), b.sort_index())
        add(5, sp, "", "rho of every regulon with pseudotime", "identical" if ok else "differs", detail)
    ok, detail = files_equal(f"{d}/regulon_quadrants_input.csv", f"{recout}/regulon_quadrants_input_fibdec_latent.csv")
    if ok is None:
        ok, detail = frames_equal(pd.read_csv(f"{d}/regulon_quadrants_input.csv"),
                                  pd.read_csv(f"{recout}/regulon_quadrants_input_fibdec_latent.csv"))
    else:
        detail += " as the table of record (sha256 6c3c99e0...)"
    add(5, "both", "", "Fig. 5E table (111 TF pairs)", "identical" if ok else "differs", detail)
    run_counts = quadrant_counts(pd.read_csv(f"{d}/regulon_quadrants_input.csv"))
    rec_counts = quadrant_counts(pd.read_csv(f"{recout}/regulon_quadrants_input_fibdec_latent.csv"))
    add(5, "both", "", "Fig. 5E counts at |rho| >= 0.2", "identical" if run_counts == rec_counts else "differs",
        ", ".join(f"{k} {v}" for k, v in run_counts.items())
        + ("" if run_counts == rec_counts else "; record: " + ", ".join(f"{k} {v}" for k, v in rec_counts.items())))


# ---------------------------------------------------------------- step 6: figures
def svg_text(path):
    s = open(path, encoding="utf-8").read()
    return [re.sub(r"\s+", " ", t).strip() for t in re.findall(r"<text[^>]*>(.*?)</text>", s, re.S)]


def step6(run, rec, recout):
    d = f"{run}/figures"
    if not os.path.isdir(d):
        add(6, "both", "", "figures", "not run")
        return
    ok, detail = files_equal(f"{d}/Fig5E_regulon_quadrants.svg",
                             f"{recout}/Fig5E_regulon_pseudotime_rho_human_vs_mouse_not_v9.svg")
    add(6, "both", "", "Fig. 5E SVG", "identical" if ok else "differs",
        (detail + " as the manuscript panel") if ok else "bytes differ from the manuscript panel")
    a, b = svg_text(f"{d}/FigS28_regulon_recurrence_cutoffs.svg"), svg_text(f"{recout}/scenic_supp.svg")
    add(6, "both", "", "Fig. S28 SVG: every text label and number", "identical" if a == b else "differs",
        f"{len(a)} labels, same order (the manuscript file was drawn with Helvetica Neue, so its bytes differ)"
        if a == b else f"{len(a)} vs {len(b)} labels; first difference "
                       f"{next(((x, y) for x, y in zip(a, b) if x != y), None)}")
    for sp in SPECIES:
        x = pd.read_csv(f"{d}/s28_inputs/{sp}_consensus_regulons_long.tsv", sep="\t")
        y = pd.read_csv(f"{rec}/scenic_seeds/{sp}/{CONS2}/consensus_regulons_long.tsv", sep="\t")
        cols = ["regulon", "target"]
        common = [c for c in y.columns if c in x.columns]
        ok, detail = frames_equal(x.sort_values(cols).reset_index(drop=True)[common],
                                  y.sort_values(cols).reset_index(drop=True)[common])
        add(6, sp, "", "Fig. S28 input: floor-2 consensus links", "identical" if ok else "differs", detail)
    a, b = json.load(open(f"{d}/s28_inputs/perrun_sizes.json")), json.load(open(f"{recout}/perrun_sizes.json"))
    for sp in SPECIES:
        ok = a[sp]["per_tf_median_size"] == b[sp]["per_tf_median_size"] and a[sp]["regulons_per_run"] == b[sp]["regulons_per_run"]
        add(6, sp, "", "Fig. S28 input: regulon sizes in single runs", "identical" if ok else "differs",
            f"regulons per run {b[sp]['regulons_per_run']}, {len(b[sp]['per_tf_median_size'])} TFs")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True, help="the run folder (config.sh OUT)")
    ap.add_argument("--record", required=True, help="the run of record's lane (config.sh RECORD)")
    ap.add_argument("--record-outputs", required=True, help="record/outputs of the pipeline folder")
    ap.add_argument("--runs", required=True, help="runs.tsv")
    ap.add_argument("--out", required=True, help="output table (TSV)")
    a = ap.parse_args()

    for sp in SPECIES:
        step1(a.run, a.record, sp)
    runs = pd.read_csv(a.runs, sep="\t")
    for sp, seed in runs.itertuples(index=False):
        step2(a.run, a.record, sp, int(seed))
    for sp, seed in runs.itertuples(index=False):
        step3(a.run, a.record, sp, int(seed))
    for sp in SPECIES:
        step4(a.run, a.record, sp)
    step5(a.run, a.record_outputs)
    step6(a.run, a.record, a.record_outputs)

    table = pd.DataFrame(ROWS, columns=["step", "species", "seed", "item", "result", "detail"])
    table.to_csv(a.out, sep="\t", index=False)
    print(f"\nwrote {a.out}\n\nsummary (comparisons per step):")
    names = {1: "pooled objects", 2: "GRNBoost2", 3: "cisTarget", 4: "consensus and AUCell", 5: "correlations",
             6: "figures"}
    for step, g in table.groupby("step"):
        c = collections.Counter(g.result)
        print(f"  {step} {names[step]:22s} " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())))
    bad = table[table.result == "differs"]
    if len(bad):
        print(f"\n{len(bad)} comparisons differ:")
        print(bad.to_string(index=False))


if __name__ == "__main__":
    main()
