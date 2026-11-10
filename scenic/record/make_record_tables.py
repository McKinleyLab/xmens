#!/usr/bin/env python
"""Writes record/SOURCES.tsv and record/step_counts.tsv for the SCENIC run of record. Read-only on every source.

SOURCES.tsv: every input, script and output of the run of record that this folder uses or compares against: what it
is, where it lives (and, for files copied here, the copy), bytes and sha256. The 100 runs per species are summarized
per folder (files, bytes), not hashed one by one.
step_counts.tsv: the numbers at every step, per species, read from the record's own files.
Run once, on the cluster, from the pipeline folder: python record/make_record_tables.py
"""
import glob
import hashlib
import json
import os
import subprocess

import numpy as np
import pandas as pd

NH = "/n/lab_storage/eddy_lab/users/nhilgert"
LANE = f"{NH}/decidualization/revisions/scenic_v2_2026_08_28"
SEEDS = f"{LANE}/scenic_seeds"
REPO = f"{NH}/xmens_final"
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAG, CONS6, CONS2 = "sketch15k", "consensus_sketch15k_weighted_soft_occ6", "consensus_sketch15k_weighted_soft"
AUCELL = CONS6 + "_pyscenic_cli"
SPECIES = ("mouse", "human")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


rows = []


def src(group, what, path, copy="", note=""):
    """One row of SOURCES.tsv; a copy in this folder must be byte-identical to its original."""
    real = os.path.realpath(path)
    digest = sha256(real)
    if copy and os.path.realpath(os.path.join(HERE, copy)) != real:
        assert sha256(os.path.join(HERE, copy)) == digest, f"{copy} differs from {path}"
    rows.append({"group": group, "what": what, "path": path, "copy_in_this_folder": copy,
                 "bytes": os.path.getsize(real), "sha256": digest, "note": note})


def folder(group, what, path, pattern, note=""):
    files = sorted(glob.glob(os.path.join(path, pattern)))
    rows.append({"group": group, "what": what, "path": path, "copy_in_this_folder": "",
                 "bytes": sum(os.path.getsize(f) for f in files), "sha256": f"{len(files)} files, not hashed",
                 "note": note})


# ---- inputs
for sp, name in (("mouse", "XMens"), ("human", "Human")):
    src("input", f"step 1 input, {sp}: fib/dec object, release 2026_09_15",
        f"{REPO}/data/deliverables/2026_09_15/Cevrim_{name}_Only_Fib_Dec_Cells.h5ad",
        note="the object SCENIC read (byte-identical to netscratch 2026_08_28_Cevrim_*_Only_Fib_Dec_Cells.h5ad)")
for sp, name in (("mouse", "XMens"), ("human", "Human")):
    src("input", f"step 5 input, {sp}: final fib/dec object (obs pseudotime)",
        f"/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects/Cevrim_{name}_Only_Fib_Dec_Cells.h5ad",
        note="sha256 as of the day this table was written; the file may be replaced")
for f in sorted(glob.glob(f"{NH}/decidualization/revisions/scenic/databases/*")):  # the lane's databases/ links here
    src("input", "cisTarget database or TF list", f)
src("input", "Ensembl 116 one-to-one pairs (step 5)", f"{REPO}/results/mi_dpt/overlap_of_record/ortholog_pairs_one2one.csv",
    copy="resources/ortholog_pairs_one2one.csv", note="written by rule mi_overlap of xmens_final")
src("input", "Ensembl 116 Compara component table (step 6)", f"{REPO}/resources/orthology/ensembl_human_mouse_components.tsv",
    note="tracked in xmens_final")

# ---- scripts of the run of record (unedited copies in record/scripts/)
for name in ("grn_only.py", "ctx_only.py", "build_consensus_aucell.py"):
    src("record script", f"{name}, as run (steps 2-4)", f"{LANE}/scripts/{name}", copy=f"record/scripts/{name}",
        note="the lane's scripts/ folder (a link to revisions/scenic/scripts)")
src("record script", "scenic_make_loom.py (steps 1-2)", f"{REPO}/workflow/scripts/scenic/scenic_make_loom.py",
    copy="record/scripts/scenic_make_loom.py", note="the record's manifests name its earlier path, workflow/scripts/")
src("record script", "24_scenic.smk, the Snakemake rules behind steps 1-4", f"{REPO}/workflow/rules/24_scenic.smk",
    copy="record/config/24_scenic.smk")
# the September rules as first committed (b9fa19c, 2026-09-14), from git
rules13 = subprocess.run(["git", "-C", REPO, "show", "b9fa19c:workflow/rules/13_scenic.smk"], check=True,
                         capture_output=True).stdout
assert hashlib.sha256(rules13).hexdigest() == sha256(f"{HERE}/record/config/13_scenic.smk"), "13_scenic.smk differs from git"
rows.append({"group": "record script", "what": "13_scenic.smk, the rules of steps 1-4 at the first commit",
             "path": f"{REPO} git b9fa19c:workflow/rules/13_scenic.smk", "copy_in_this_folder": "record/config/13_scenic.smk",
             "bytes": len(rules13), "sha256": hashlib.sha256(rules13).hexdigest(), "note": "committed 2026-09-14"})
R5 = f"{NH}/tmp/2026_10_30_5e_samap_s30fg/item1_5e"
for name, what, orig in (("scenic_quadrant_rho_fibdec_latent.py", "the Fig. 5E table of record, 30 Oct 2026",
                          f"{R5}/new/scenic_quadrant_rho_fibdec_latent.py"),
                         ("plot_regulon_quadrants.py", "the manuscript's Fig. 5E, 30 Oct 2026", f"{R5}/plot_regulon_quadrants.py")):
    src("record script", f"{name}: {what}", orig, copy=f"record/scripts/{name}")
for name, what in (("consensus_floor6_job.sh", "the floor-6 consensus job, 14 Oct 2026 (Slurm 52987211)"),
                   ("cliauc_run.sh", "the AUCell of record, 26 Oct 2026"),
                   ("regulon_recurrence_cutoffs_fig.py", "the manuscript's Fig. S28, 18 Oct - 2 Nov 2026"),
                   ("perrun_regulon_sizes.py", "the per-run regulon sizes of Fig. S28, 17 Sep 2026")):
    src("record script", f"{name}: {what}", f"{HERE}/record/scripts/{name}", copy=f"record/scripts/{name}",
        note="kept in the Claude Science project of the analysis; this folder holds the only cluster copy")

# ---- outputs of the run of record
for sp in SPECIES:
    for ext in ("h5ad", "loom", "manifest.json"):
        src("record output", f"step 1, {sp}: {sp}_dsc_pooled.{ext}", f"{SEEDS}/{sp}/{sp}_dsc_pooled.{ext}")
    folder("record output", f"step 2, {sp}: the runs' sketches", f"{SEEDS}/{sp}/looms_{TAG}", "*.loom")
    folder("record output", f"step 2, {sp}: the runs' adjacencies", f"{SEEDS}/{sp}/adjacencies_{TAG}", "*.tsv")
    folder("record output", f"step 3, {sp}: the runs' regulons", f"{SEEDS}/{sp}/regulons_{TAG}", "*.pkl")
    for cons, floor in ((CONS6, 6), (CONS2, 2)):
        for name in ("recurrence_report.tsv", "consensus_regulons.pkl", "auc_mtx.csv"):
            src("record output", f"step 4, {sp}: floor {floor} {name}", f"{SEEDS}/{sp}/{cons}/{name}")
    src("record output", f"step 6 input, {sp}: floor-2 links (10 Sep 2026)", f"{SEEDS}/{sp}/{CONS2}/consensus_regulons_long.tsv")
    for name in ("consensus_regulons.yaml", "auc_mtx_weighted.csv", "auc_mtx_unweighted.csv"):
        src("record output", f"step 4, {sp}: AUCell of record {name}", f"{SEEDS}/{sp}/{AUCELL}/{name}")
for name, what, orig in (
        ("regulon_quadrants_input_fibdec_latent.csv", "step 5: the Fig. 5E table of record",
         f"{REPO}/resources/figures/panels_programs/regulon_quadrants_input_fibdec_latent.csv"),
        ("regulon_rho_pseudotime_mouse_fibdec_latent.tsv", "step 5: rho of every mouse regulon",
         f"{NH}/tmp/2026_10_30_5e_samap_s30fg/item1_5e/new/regulon_rho_pseudotime_mouse_fibdec_latent.tsv"),
        ("regulon_rho_pseudotime_human_fibdec_latent.tsv", "step 5: rho of every human regulon",
         f"{NH}/tmp/2026_10_30_5e_samap_s30fg/item1_5e/new/regulon_rho_pseudotime_human_fibdec_latent.tsv"),
        ("Fig5E_regulon_pseudotime_rho_human_vs_mouse_not_v9.svg", "step 6: the manuscript's Fig. 5E",
         f"{NH}/xmens_figures/final_figures/Fig5_human_mouse_shared_programs/Fig5E_regulon_pseudotime_rho_human_vs_mouse_not_v9.svg"),
        ("scenic_supp.svg", "step 6: the manuscript's Fig. S28",
         f"{NH}/xmens_figures/final_figures/S28_scenic_supp/scenic_supp.svg"),
        ("perrun_sizes.json", "step 6 input: regulon sizes in single runs (2 Nov 2026)",
         f"{REPO}/results/figures/panels/scenic_params/perrun_sizes.json")):
    copy = f"{HERE}/record/outputs/{name}"
    assert sha256(copy) == sha256(orig), f"copy of {orig} differs from the original"
    src("record output", what, orig, copy=f"record/outputs/{name}")
pd.DataFrame(rows).to_csv(f"{HERE}/record/SOURCES.tsv", sep="\t", index=False)
print(f"SOURCES.tsv: {len(rows)} rows", flush=True)

# ---- step counts of the run of record
counts = []


def count(step, quantity, values, source):
    counts.append({"step": step, "quantity": quantity, "mouse": values["mouse"], "human": values["human"], "source": source})


man = {sp: json.load(open(f"{SEEDS}/{sp}/{sp}_dsc_pooled.manifest.json")) for sp in SPECIES}
for key, q in (("n_cells_full", "nuclei"), ("n_genes_full", "genes"), ("gene_filter_min_cells", "nuclei a gene must be detected in (1%)"),
               ("n_genes_grn", "genes in the pooled loom"), ("n_tfs_in_list", "TFs in the TF list"),
               ("n_tfs_retained", "TFs among the loom's genes")):
    count(1, q, {sp: man[sp][key] for sp in SPECIES}, "<species>_dsc_pooled.manifest.json")

links, sketch_n = {}, {}
for sp in SPECIES:
    n = []
    for f in sorted(glob.glob(f"{SEEDS}/{sp}/adjacencies_{TAG}/adjacencies_seed*.tsv")):
        with open(f, "rb") as fh:
            n.append(sum(1 for _ in fh) - 1)
    links[sp] = n
    sketch_n[sp] = sorted({json.load(open(f))["n_cells_grn"] for f in glob.glob(f"{SEEDS}/{sp}/looms_{TAG}/*.manifest.json")})
count(2, "runs (seeds 1000-1099)", {sp: len(links[sp]) for sp in SPECIES}, "adjacencies_sketch15k/")
count(2, "nuclei per sketch", {sp: ",".join(map(str, sketch_n[sp])) for sp in SPECIES}, "looms_sketch15k/*.manifest.json")
count(2, "links per run (min / median / max)",
      {sp: f"{min(links[sp]):,} / {int(np.median(links[sp])):,} / {max(links[sp]):,}" for sp in SPECIES}, "adjacencies_sketch15k/")
per = json.load(open(f"{HERE}/record/outputs/perrun_sizes.json"))
count(3, "regulons per run (min / median / max)",
      {sp: f"{per[sp]['regulons_per_run']['min']} / {per[sp]['regulons_per_run']['median']:g} / {per[sp]['regulons_per_run']['max']}"
       for sp in SPECIES}, "record/outputs/perrun_sizes.json")
count(3, "TFs with a regulon in at least one run", {sp: per[sp]["n_distinct_tfs_any_run"] for sp in SPECIES},
      "record/outputs/perrun_sizes.json")
for cons, floor in ((CONS6, 6), (CONS2, 2)):
    rep = {sp: pd.read_csv(f"{SEEDS}/{sp}/{cons}/recurrence_report.tsv", sep="\t") for sp in SPECIES}
    count(4, f"floor {floor}: regulons kept (in >= 80 of 100 runs, >= 10 targets)", {sp: int(rep[sp].kept.sum()) for sp in SPECIES},
          f"{cons}/recurrence_report.tsv")
    count(4, f"floor {floor}: median targets per kept regulon",
          {sp: int(np.median(rep[sp][rep[sp].kept].n_targets_kept)) for sp in SPECIES}, f"{cons}/recurrence_report.tsv")
shape = {}
for sp in SPECIES:
    with open(f"{SEEDS}/{sp}/{AUCELL}/auc_mtx_weighted.csv") as fh:
        n_reg = len(fh.readline().rstrip("\n").split(",")) - 1
        shape[sp] = f"{sum(1 for _ in fh):,} x {n_reg}"
count(4, "activity of record (nuclei x regulons)", shape, f"{AUCELL}/auc_mtx_weighted.csv")
rho = {sp: pd.read_csv(f"{HERE}/record/outputs/regulon_rho_pseudotime_{sp}_fibdec_latent.tsv", sep="\t") for sp in SPECIES}
count(5, "regulons with a rho", {sp: len(rho[sp]) for sp in SPECIES}, "record/outputs/regulon_rho_pseudotime_*_fibdec_latent.tsv")
t = pd.read_csv(f"{HERE}/record/outputs/regulon_quadrants_input_fibdec_latent.csv")
h, m = t.rho_human.to_numpy(), t.rho_mouse.to_numpy()
both = (np.abs(h) >= 0.2) & (np.abs(m) >= 0.2)
fig = {"TF pairs (Ensembl one-to-one, a regulon in both species)": len(t),
       "pairs with |rho| >= 0.2 in both species": int(both.sum()),
       "up in both": int((both & (h > 0) & (m > 0)).sum()), "down in both": int((both & (h < 0) & (m < 0)).sum()),
       "up in human, down in mouse": int((both & (h > 0) & (m < 0)).sum()),
       "down in human, up in mouse": int((both & (h < 0) & (m > 0)).sum())}
for q, v in fig.items():
    counts.append({"step": 5, "quantity": q, "mouse": v, "human": v, "source": "record/outputs/regulon_quadrants_input_fibdec_latent.csv"})
pd.DataFrame(counts).to_csv(f"{HERE}/record/step_counts.tsv", sep="\t", index=False)
print(pd.DataFrame(counts).to_string(index=False))
