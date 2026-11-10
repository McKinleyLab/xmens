#!/usr/bin/env python
"""Inputs for every pseudotime heatmap, read from the final objects and the MI of record.

Objects: final_objects/Cevrim_{XMens,Human}_Only_Fib_Dec_Cells.h5ad (2 Nov 2026). obs['pseudotime'] is the
rank-normalized pseudotime of record: mouse = the 2026_10_22 harmonized pseudotime (pt_harmonized_rank), human = the
30 Oct 2026 refit (data/interim/human_pt_2026_10_30_run4_a40/cells.parquet, pt_rank). Both are checked here against
those tables. Nuclei with NaN pseudotime (mouse decAp, by design) are dropped. Raw counts: layers['counts'].

Gene sets
  * Yang et al. (2023) 28-gene panel (mouse only).
  * MI heatmaps: Fig. S16 top-decile sets = genes whose MI exceeds the 90th percentile of MI among each species'
    significant genes (q < 0.05), restricted to strict one-to-one orthologs significant in both species
    (procedure of prep_fullheat_args.py). Checked against overlap_of_record/methods_overlap.json.
  * Cross-species heatmap: the crossspecies_gene_selection.py procedure (pool = top 10% by MI of both species over the
    one-to-one orthologs measured in both MI arms; 50-bin log1p-CP10K z-profiles; r > 0.8; 14 early + 14 late by r;
    rows ordered by mean centroid), run on the final objects' pseudotime. Checked against the pipeline's selection
    (overlap_of_record/crossspecies_heatmap_selection/conserved_profiles_selected_genes.csv).
Dispersion: Monocle 2 disp_calc_helper_NB (size factors = total counts / geometric mean; genes detected in more than
one nucleus; moments estimate with the sample variance, n - 1; negatives set to 0) on the nuclei with pseudotime.

usage: export_heatmap_inputs_final_objects.py --mouse-h5ad H5AD --human-h5ad H5AD --mi-dir DIR --pairs CSV
           --pt-record-mouse PARQUET --pt-record-human PARQUET --record-overlap JSON --record-selection CSV
           --previous-selection CSV --out DIR
"""
import argparse, hashlib, json, os, time
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy.stats import hypergeom, spearmanr

ap = argparse.ArgumentParser()
for name in ("--mouse-h5ad", "--human-h5ad", "--mi-dir", "--pairs", "--pt-record-mouse", "--pt-record-human",
             "--record-overlap", "--record-selection", "--previous-selection", "--out"):
    ap.add_argument(name, required=True)
ap.add_argument("--mouse-arm", default="pt_harmonized_mouse_nall"); ap.add_argument("--human-arm", default="pt_harmonized_human_nall")
ap.add_argument("--topfrac", type=float, default=0.10); ap.add_argument("--rmin", type=float, default=0.8)
ap.add_argument("--n-early", type=int, default=14); ap.add_argument("--n-late", type=int, default=14)
a = ap.parse_args()
for sub in ("", "yang", "mi", "crossspecies"):
    os.makedirs(os.path.join(a.out, sub), exist_ok=True)
SPECIES = ("mouse", "human")
H5 = {"mouse": a.mouse_h5ad, "human": a.human_h5ad}
ARMS = {"mouse": a.mouse_arm, "human": a.human_arm}
PTREC = {"mouse": a.pt_record_mouse, "human": a.pt_record_human}
YANG = ["FN1", "VIM", "COL3A1", "COL1A2", "COL1A1", "IGFBP5", "ZEB1", "SFRP4", "IGF1", "TCF4", "CDH11", "FSCN1", "WNT5A",
        "F2R", "ASS1", "UBE2C", "TOP2A", "MKI67", "HMGB1", "JUN", "BMP2", "LAMA5", "ANGPT4", "ICAM1", "ADAM19", "TREM2", "CDH5", "GJA1"]
DAY = {"uninduced": 0.0, "1.5_dpi_early": 1.5, "1.5_dpi_late": 1.5, "2.5_dpi": 2.5, "3.5_dpi": 3.5, "4.5_dpi_bleeding": 4.5}
# day encoding stated in the Methods (Pseudotime section)
DOC_DAYS = {"mouse": {"uninduced": 0.0, "1.5_dpi_early": 1.0, "1.5_dpi_late": 1.5, "2.5_dpi": 2.0, "3.5_dpi": 3.0, "4.5_dpi_bleeding": 4.0},
            "human": {"secretory": 0.0, "menstrual": 1.0}}
T0 = time.time()
def log(*x): print(f"[{time.time() - T0:7.1f}s]", *x, flush=True)
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()
def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (str, bool, int, float)) or x is None:
        return x
    if isinstance(x, (np.integer, np.floating, np.bool_)):
        return x.item()
    if isinstance(x, np.ndarray) and x.size <= 20:
        return x.tolist()
    return f"<{type(x).__name__}>"
CHECKS = dict(inputs={k: dict(path=v, sha256=sha256(v)) for k, v in dict(
    mouse_h5ad=a.mouse_h5ad, human_h5ad=a.human_h5ad, pairs=a.pairs, pt_record_mouse=a.pt_record_mouse,
    pt_record_human=a.pt_record_human, record_overlap=a.record_overlap, record_selection=a.record_selection,
    previous_selection=a.previous_selection, mi_mouse=f"{a.mi_dir}/{a.mouse_arm}/mi_all.csv",
    mi_human=f"{a.mi_dir}/{a.human_arm}/mi_all.csv").items()})
log("hashed inputs")

# ---- 1. objects: nuclei with pseudotime, raw counts, labels; pseudotime checked against the record tables -------
D = {}
for s in SPECIES:
    obj = ad.read_h5ad(H5[s])
    log(s, obj.shape)
    pt_all = pd.to_numeric(obj.obs["pseudotime"], errors="coerce").to_numpy(float)
    keep = np.isfinite(pt_all)
    C = obj.layers["counts"]
    C = (C if sp.issparse(C) else sp.csr_matrix(C)).tocsr()[keep].astype(np.float64)
    obs = obj.obs.loc[keep, ["cell_subtype", "timepoint", "sample"]].astype(str)
    obs.index = pd.Index(obj.obs_names.astype(str))[keep]
    info = dict(n_obs=int(obj.n_obs), n_vars=int(obj.n_vars), n_with_pseudotime=int(keep.sum()),
                nan_by_subtype=obj.obs.loc[~keep, "cell_subtype"].astype(str).value_counts().to_dict(),
                counts_integer=bool(np.allclose(C.data, np.round(C.data))),
                n_unique_pseudotime=int(len(np.unique(pt_all[keep]))),
                pseudotime_min_max=[float(pt_all[keep].min()), float(pt_all[keep].max())],
                uns_pseudotime=jsonable(obj.uns.get("pseudotime", {})),
                subtypes=obs.cell_subtype.value_counts().to_dict(), timepoints=obs.timepoint.value_counts().to_dict())
    var = pd.Index(obj.var_names.astype(str))
    del obj
    rec = pd.read_parquet(PTREC[s])
    rec.index = rec.index.astype(str)
    rec_pt = pd.to_numeric(rec["pt_rank"], errors="coerce")
    rec_finite = set(rec_pt.index[rec_pt.notna()])
    aligned = rec_pt.reindex(obs.index).to_numpy(float)
    info["record_check"] = dict(table=PTREC[s], n_record=int(len(rec)), n_record_finite=len(rec_finite),
                                same_nuclei=rec_finite == set(obs.index),
                                max_abs_diff=float(np.nanmax(np.abs(aligned - pt_all[keep]))),
                                n_missing_in_record=int(np.isnan(aligned).sum()))
    log(s, "record check", info["record_check"])
    # Methods check: pseudotime vs timepoint with the Methods' day encoding (mouse rho = 0.62, human rho = 0.32)
    days = obs.timepoint.map(DOC_DAYS[s]).to_numpy(float)
    info["spearman_pt_vs_timepoint_methods_encoding"] = dict(
        rho=round(float(spearmanr(pt_all[keep][np.isfinite(days)], days[np.isfinite(days)]).correlation), 4),
        n=int(np.isfinite(days).sum()), unmapped=sorted(set(obs.timepoint[~np.isfinite(days)])))
    total = np.asarray(C.sum(axis=1)).ravel()
    D[s] = dict(C=C, pt=pt_all[keep], obs=obs, var=var, total=total)
    CHECKS[s] = info

def percell(s, genes):
    d = D[s]
    cols = d["var"].get_indexer(pd.Index(genes))
    assert (cols >= 0).all(), f"{s}: genes missing from object: {[g for g, c in zip(genes, cols) if c < 0]}"
    G = np.round(d["C"][:, cols].toarray()).astype(np.int32)
    df = pd.DataFrame(G, columns=[f"cnt_{g}" for g in genes], index=d["obs"].index)
    df.insert(0, "total_counts", np.round(d["total"]).astype(np.int64)); df.insert(1, "pt_rank_norm", d["pt"])
    df.index.name = "cell"
    return df

# ---- 2. dispersion tables and cell tables ---------------------------------------------------------------------------
for s in SPECIES:
    d = D[s]
    R = d["C"].copy(); R.data = np.round(R.data)
    ndet = np.asarray((R > 0.1).sum(axis=0)).ravel(); nz = ndet > 1
    sf = d["total"] / np.exp(np.mean(np.log(d["total"])))
    Xn = (sp.diags(1.0 / sf) @ R[:, nz]).tocsc(); n = Xn.shape[0]
    mu = np.asarray(Xn.mean(axis=0)).ravel()
    ex2 = np.asarray(Xn.multiply(Xn).mean(axis=0)).ravel()
    var = (ex2 - mu ** 2) * n / (n - 1)
    disp = (var - float(np.mean(1.0 / sf)) * mu) / mu ** 2; disp[disp < 0] = 0
    pd.DataFrame({"gene": d["var"][nz], "mu": mu, "disp": disp, "n_detected": ndet[nz]}).to_csv(
        f"{a.out}/disp_table_{s}.tsv.gz", sep="\t", index=False, compression="gzip")
    CHECKS[s]["disp_n_genes"] = int(nz.sum())
    del R, Xn
    cells = pd.DataFrame({"cell": d["obs"].index, "pt_rank_norm": d["pt"], "label": d["obs"].cell_subtype.to_numpy(),
                          "timepoint": d["obs"].timepoint.to_numpy()})
    cells.to_csv(f"{a.out}/cells_{s}.tsv.gz", sep="\t", index=False, compression="gzip")
    log(s, "disp table", int(nz.sum()), "genes; cells table", cells.shape)

# ---- 3. Yang panel per-nucleus table (format of export_percell_pt_harmonized.py; 'pt' = the rank pseudotime, the
#         final objects carry no raw DPT column and the smoother only reads pt_rank_norm) -----------------------------
d = D["mouse"]
up = {g.upper(): g for g in d["var"]}
present = [g for g in YANG if g.upper() in up]
CHECKS["yang"] = dict(genes_missing=[g for g in YANG if g not in present])
P = percell("mouse", [up[g.upper()] for g in present])
y = pd.DataFrame({"cell": P.index, "pt": d["pt"], "total_counts": P.total_counts.to_numpy(),
                  "timepoint": d["obs"].timepoint.to_numpy(), "sample": d["obs"]["sample"].to_numpy(),
                  "cell_subtype": d["obs"].cell_subtype.to_numpy(), "pt_rank_norm": d["pt"]})
for g in present:
    cnt = P[f"cnt_{up[g.upper()]}"].to_numpy()
    y[f"cnt_{g}"] = cnt
    y[f"log_{g}"] = np.log1p(cnt / np.maximum(y.total_counts.to_numpy(), 1) * 1e4).round(5)
y.to_csv(f"{a.out}/yang/percell_genes_final.tsv.gz", sep="\t", index=False, compression="gzip")
day = y.timepoint.map(DAY)
CHECKS["yang"].update(shape=list(y.shape), n_induced=int((y.timepoint != "uninduced").sum()),
                      timepoints=y.timepoint.value_counts().to_dict(),
                      spearman_pt_vs_day=round(float(spearmanr(y.pt_rank_norm[day.notna()], day[day.notna()]).correlation), 4),
                      unmapped_timepoints=sorted(set(y.timepoint) - set(DAY)))
log("yang table", y.shape)

# ---- 4. MI top-decile sets (Fig. S16 definition) and their per-nucleus tables ----------------------------------------
mi = {s: pd.read_csv(f"{a.mi_dir}/{ARMS[s]}/mi_all.csv").drop_duplicates("gene").set_index("gene") for s in SPECIES}
pairs_all = pd.read_csv(a.pairs)
sig = {s: mi[s][mi[s].q < 0.05] for s in SPECIES}
bg = pairs_all[pairs_all.mouse.isin(sig["mouse"].index) & pairs_all.human.isin(sig["human"].index)].copy()
bg["mi_m"] = sig["mouse"].mi.reindex(bg.mouse).to_numpy(); bg["mi_h"] = sig["human"].mi.reindex(bg.human).to_numpy()
thr = {s: float(np.percentile(sig[s].mi, 90)) for s in SPECIES}
tm, th = bg.mi_m > thr["mouse"], bg.mi_h > thr["human"]
N, K, n_h, k = len(bg), int(tm.sum()), int(th.sum()), int((tm & th).sum()); e = K * n_h / N
mi_summary = dict(arms=ARMS, tested={s: len(mi[s]) for s in SPECIES}, significant={s: len(sig[s]) for s in SPECIES},
                  thr_nats=thr, top_genomewide={s: int((sig[s].mi > thr[s]).sum()) for s in SPECIES},
                  N=N, K_mouse=K, n_human=n_h, shared=k, expected=round(e, 3), fold=round(k / e, 4),
                  log10_p=round(float(hypergeom.logsf(k - 1, N, K, n_h) / np.log(10)), 3))
rec = json.load(open(a.record_overlap))
mi_summary["record_overlap"] = {key: rec.get(key) for key in ("N", "K", "n", "k", "expected", "fold", "log10_p")}
mi_summary["matches_record"] = bool((N, K, n_h, k) == (rec["N"], rec["K"], rec["n"], rec["k"]))
sets = {"mouse": bg[tm].assign(gene=bg.mouse, mi=bg.mi_m, partner=bg.human, shared=(tm & th)[tm]),
        "human": bg[th].assign(gene=bg.human, mi=bg.mi_h, partner=bg.mouse, shared=(tm & th)[th])}
for s in SPECIES:
    t = sets[s][["gene", "partner", "mi", "shared"]].sort_values("mi", ascending=False)
    t.to_csv(f"{a.out}/mi/topdecile_{s}.csv", index=False)
    P = percell(s, list(t.gene))
    P.to_csv(f"{a.out}/mi/percell_{s}_topdecile.tsv.gz", sep="\t", compression="gzip")
    mi_summary[f"percell_{s}"] = list(P.shape)
    log(s, "MI top-decile table", P.shape)
json.dump(mi_summary, open(f"{a.out}/mi/summary.json", "w"), indent=1)
CHECKS["mi"] = mi_summary
log("MI sets", json.dumps({kk: mi_summary[kk] for kk in ("N", "K_mouse", "n_human", "shared", "matches_record")}))

# ---- 5. cross-species selection (crossspecies_gene_selection.py procedure) on the final pseudotime --------------------
miv = {s: mi[s]["mi"] for s in SPECIES}
pairs = pairs_all[["human", "mouse"]]
pairs = pairs[pairs.human.isin(miv["human"].index) & pairs.mouse.isin(miv["mouse"].index)].reset_index(drop=True)
mm, hh = miv["mouse"].reindex(pairs.mouse).to_numpy(), miv["human"].reindex(pairs.human).to_numpy()
ok = np.isfinite(mm) & np.isfinite(hh); pairs, mm, hh = pairs[ok].reset_index(drop=True), mm[ok], hh[ok]
ktop = int(round(a.topfrac * len(pairs)))
in_top = (pd.Series(mm).rank(ascending=False).to_numpy() <= ktop) & (pd.Series(hh).rank(ascending=False).to_numpy() <= ktop)
pool = pairs[in_top].assign(mi_mouse=mm[in_top], mi_human=hh[in_top]).reset_index(drop=True)
# Methods check: the pool is "a stricter subset of the 373 shared top genes"
shared_pairs = set(zip(bg.mouse[tm & th], bg.human[tm & th]))
pool_in_shared = int(sum((m_, h_) in shared_pairs for m_, h_ in zip(pool.mouse, pool.human)))
PM, PH = percell("mouse", list(pool.mouse)), percell("human", list(pool.human))
PM.to_csv(f"{a.out}/crossspecies/percell_mouse_pool.tsv.gz", sep="\t"); PH.to_csv(f"{a.out}/crossspecies/percell_human_pool.tsv.gz", sep="\t")
def profiles(P, genes, nbin):
    X = np.log1p(P[[f"cnt_{g}" for g in genes]].to_numpy(float) / P["total_counts"].to_numpy(float)[:, None] * 1e4)
    q = pd.qcut(P["pt_rank_norm"].rank(method="first"), nbin, labels=False).to_numpy()
    Mx = np.vstack([X[q == b].mean(axis=0) for b in range(nbin)]).T
    sd = Mx.std(axis=1, ddof=0); sd[sd == 0] = 1
    return Mx, (Mx - Mx.mean(axis=1, keepdims=True)) / sd[:, None]
def centroid(Z):
    W = np.clip(Z, 0, None); pos = np.arange(Z.shape[1]) / (Z.shape[1] - 1)
    return (W * pos).sum(1) / W.sum(1)
corr = lambda A, B: np.array([np.corrcoef(A[i], B[i])[0, 1] for i in range(len(A))])
Mh50, Zh50 = profiles(PH, list(pool.human), 50); Mm50, Zm50 = profiles(PM, list(pool.mouse), 50)
_, Zh20 = profiles(PH, list(pool.human), 20); _, Zm20 = profiles(PM, list(pool.mouse), 20)
T = pool.copy()
T["r50"], T["r20"] = corr(Zh50, Zm50), corr(Zh20, Zm20)
T["centroid_human"], T["centroid_mouse"] = centroid(Zh50), centroid(Zm50)
T["pt_centroid"] = (T.centroid_human + T.centroid_mouse) / 2
T["peak_h"], T["peak_m"] = Zh50.argmax(1), Zm50.argmax(1)
T["range_h"], T["range_m"] = Mh50.max(1) - Mh50.min(1), Mm50.max(1) - Mm50.min(1)
T["phase"] = np.where(T.pt_centroid < 0.4, "early", np.where(T.pt_centroid > 0.6, "late", "mid"))
T.sort_values("r50", ascending=False).to_csv(f"{a.out}/crossspecies/conserved_profiles_pool_concordance.csv", index=False)
hi = T[T.r50 > a.rmin]
pick = pd.concat([hi[hi.phase == "early"].nlargest(a.n_early, "r50"), hi[hi.phase == "late"].nlargest(a.n_late, "r50")]).sort_values("pt_centroid")
S = pick.rename(columns={"r50": "profile_r_50bins", "r20": "profile_r_20bins"})[
    ["human", "mouse", "profile_r_50bins", "profile_r_20bins", "phase", "pt_centroid", "peak_h", "peak_m"]]
S.to_csv(f"{a.out}/crossspecies/conserved_profiles_selected_genes.csv", index=False)
percell("mouse", list(S.mouse)).to_csv(f"{a.out}/crossspecies/percell_mouse_sel28.tsv.gz", sep="\t")
percell("human", list(S.human)).to_csv(f"{a.out}/crossspecies/percell_human_sel28.tsv.gz", sep="\t")
recsel = pd.read_csv(a.record_selection); prev = pd.read_csv(a.previous_selection)
m = S.merge(recsel, on=["human", "mouse"], suffixes=("", "_record"))
sel_summary = dict(mi_arms=ARMS, universe=int(len(pairs)), k_top=ktop, pool=int(len(pool)),
                   pool_in_shared_topdecile=pool_in_shared, shared_topdecile=len(shared_pairs), n_r_gt=int(len(hi)), rmin=a.rmin,
                   phase_counts_r_gt=hi.phase.value_counts().to_dict(), median_r50_pool=round(float(T.r50.median()), 3),
                   n_selected=int(len(S)), r50_range_selected=[round(float(S.profile_r_50bins.min()), 4), round(float(S.profile_r_50bins.max()), 4)],
                   same_genes_as_record=set(S.human) == set(recsel.human),
                   same_order_as_record=list(S.human) == list(recsel.human),
                   max_abs_diff_vs_record=dict(r50=float((m.profile_r_50bins - m.profile_r_50bins_record).abs().max()),
                                               r20=float((m.profile_r_20bins - m.profile_r_20bins_record).abs().max()),
                                               centroid=float((m.pt_centroid - m.pt_centroid_record).abs().max())),
                   previous_selection=a.previous_selection, prev_columns=list(prev.columns),
                   overlap_with_previous=int(len(set(S.human) & set(prev[prev.columns[0]].astype(str)))),
                   order=list(S.human))
json.dump(sel_summary, open(f"{a.out}/crossspecies/selection_summary.json", "w"), indent=1)
CHECKS["crossspecies"] = sel_summary
log("selection", json.dumps({kk: sel_summary[kk] for kk in ("pool", "n_r_gt", "same_genes_as_record", "same_order_as_record", "max_abs_diff_vs_record")}))
json.dump(CHECKS, open(f"{a.out}/checks.json", "w"), indent=1, default=str)
log("done")
