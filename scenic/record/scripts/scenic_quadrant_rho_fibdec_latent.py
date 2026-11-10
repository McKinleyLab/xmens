#!/usr/bin/env python
"""Step 1 of the SCENIC quadrant figure: regulon activity vs the pseudotime of record, per regulon, both species.

_fibdec_latent variant (2026-10-30). Identical to scenic_quadrant_rho.py (sha256 19996b6d...) except for the HUMAN
pseudotime: pt_harmonized_rank is replaced by pt_rank of the side table recomputed on the human stromal scVI latent
refit on the 47,533 fibroblast/decidual nuclei only (HUMAN_PT_TABLE below), joined on obs_names. The side table must
hold exactly the barcodes of the human object and of the v9 release object (X/data/deliverables/2026_10_30); nothing
is re-ranked. Mouse is read exactly as in the run of record. Outputs carry the suffix _fibdec_latent.

For every floor-6 consensus regulon, Spearman rho between its per-nucleus activity and obs['pt_harmonized_rank'],
then the 5-column input table of the cross-species quadrant plot (step 2, fig_regulon_quadrants.py).

Inputs (cluster paths)
  activity   {LANE}/scenic_seeds/{sp}/consensus_sketch15k_weighted_soft_occ6_pyscenic_cli/auc_mtx_weighted.csv
             floor-6 consensus regulons (TF in >= 80 of 100 runs, targets in >= 6 runs, occurrence-weighted), scored
             with pySCENIC 0.12.1 `pyscenic aucell <filtered pooled loom> consensus_regulons.yaml --weights --seed 0`
  pseudotime pseudotime of record (2026-10-27): obs['pt_harmonized_rank'] of the fibroblast/decidual objects in
             /n/netscratch/eddy_lab/Everyone/nhilgert/2026_10_22_flagged_harmonized/. Mouse decidual_apoptotic nuclei
             are NaN by design and are dropped, not filled. Nuclei are matched to the activity matrix by obs_names; the
             script stops unless every nucleus with pseudotime has an activity score.
  orthologs  ortholog_pairs_one2one.csv: Ensembl one-to-one human-mouse symbol pairs (copy in this folder)
  labels     the 20 canonical decidualization TFs named on the plot (LABELLED; or edit the csv `label` column)

Sanity checks from the pseudotime note are recomputed and written to the provenance file (expected: 42,290 mouse and
47,533 human nuclei with pseudotime; mouse decidual medians boundary 0.543 < proliferative 0.721 < mature 0.947).
Spearman rho = Pearson correlation of average ranks; activity cast to float32 before ranking.
Outputs (working directory): regulon_rho_pseudotime_{mouse,human}.tsv, regulon_quadrants_input.csv,
scenic_quadrant_rho.provenance.json.
"""
import hashlib, json, os, sys, datetime
import numpy as np, pandas as pd, h5py
from scipy.stats import rankdata

LANE = "/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/scenic_v2_2026_08_28"
PTDIR = "/n/netscratch/eddy_lab/Everyone/nhilgert/2026_10_22_flagged_harmonized"
OBJ = {"mouse": f"{PTDIR}/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad", "human": f"{PTDIR}/Cevrim_Human_Only_Fib_Dec_Cells.h5ad"}
PTCOL = "pt_harmonized_rank"
PT_SCRIPT_SHA = "d4b618ae2b037a833e927948c4095cce19f6028b8123dfda68e60bb16553ffa1"
EXPECTED_N = {"mouse": 42290, "human": 47533}
AUC = "consensus_sketch15k_weighted_soft_occ6_pyscenic_cli/auc_mtx_weighted.csv"
ORTHO = sys.argv[1] if len(sys.argv) > 1 else "ortholog_pairs_one2one.csv"
# Human pseudotime of 2026-10-30 (fib/dec-only latent), keyed by obs_names; pt_rank = new pt_harmonized_rank.
HUMAN_PT_TABLE = ("/n/netscratch/eddy_lab/Everyone/nhilgert/2026_10_30_human_fibdec_latent/run4_a40/fibdec/"
                  "pt_human_pca_released/cells.parquet")
HUMAN_PT_COL = "pt_rank"
V9_HUMAN = ("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/deliverables/2026_10_30/"
            "Cevrim_Human_Only_Fib_Dec_Cells.h5ad")
SUFFIX = "_fibdec_latent"
LABELLED = ["CEBPB", "CEBPD", "CREB1", "ETS1", "FOS", "FOSL2", "FOXO1", "FOXO3", "GATA2", "HOXA10", "IRF1", "JUN",
            "JUNB", "KLF4", "KLF6", "NFKB1", "NR2F2", "STAT1", "TFAP2C", "WT1"]


def sha256(p, chunk=1 << 24):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def strs(ds):
    try:
        return np.asarray(ds.asstr()[:], dtype=object)
    except Exception:
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in ds[:]], dtype=object)


def obs_col(obs, name):
    node = obs[name]
    if isinstance(node, h5py.Group):                       # categorical (anndata >= 0.8)
        cats, codes = strs(node["categories"]), node["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], None)
    v = node[:]
    return v.astype(float) if v.dtype.kind in "fiu" else strs(node)


def read_obs(path):
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        idx = strs(obs[obs.attrs["_index"]])
        pt = obs_col(obs, PTCOL).astype(float)
        cols = {}
        for c in obs.keys():
            if c == obs.attrs["_index"] or c == PTCOL:
                continue
            try:
                v = obs_col(obs, c)
            except Exception:
                continue
            if v.dtype == object and np.isin(["decidual_mature", "decidual_boundary"], v).all():
                cols["subtype"] = (c, v)
            if c.lower() in ("day", "days", "day_num", "timepoint_day") and v.dtype.kind in "fiu":
                cols["day"] = (c, v)
    return pd.Index(idx), pt, cols


def spearman(A, pt):
    """Spearman rho of each column of A (cells x regulons) with pt."""
    pr = rankdata(pt); pr -= pr.mean()
    RA = rankdata(np.asarray(A, dtype=np.float32), axis=0); RA -= RA.mean(axis=0)
    return (pr @ RA) / (np.sqrt((RA * RA).sum(axis=0)) * np.sqrt((pr * pr).sum()))


prov = {"written": datetime.datetime.now().isoformat(timespec="seconds"), "pseudotime_column": {"mouse": PTCOL, "human": f"{HUMAN_PT_COL} of HUMAN_PT_TABLE (new pt_harmonized_rank)"}, "inputs": {}, "species": {}}
man = os.path.join(PTDIR, "MANIFEST.json")
prov["pt_script_sha256_in_manifest"] = (PT_SCRIPT_SHA in open(man).read()) if os.path.exists(man) else "MANIFEST.json not found"
rho = {}
for sp in ("mouse", "human"):
    idx, pt, cols = read_obs(OBJ[sp])
    if sp == "human":
        side = pd.read_parquet(HUMAN_PT_TABLE)
        with h5py.File(V9_HUMAN, "r") as f:
            v9_idx = pd.Index(strs(f["obs"][f["obs"].attrs["_index"]]))
        for name, other in ((OBJ[sp], idx), (V9_HUMAN, v9_idx)):
            assert side.index.is_unique and len(side) == len(other) and side.index.isin(other).all(), \
                f"{HUMAN_PT_TABLE} and {name} hold different barcodes"
        pt_old = pt
        pt = side.loc[idx, HUMAN_PT_COL].to_numpy(dtype=float)
        prov["human_pt_side_table"] = {
            "path": HUMAN_PT_TABLE, "sha256": sha256(HUMAN_PT_TABLE), "column": HUMAN_PT_COL,
            "run_record": os.path.join(os.path.dirname(HUMAN_PT_TABLE), "provenance.json"),
            "run_record_sha256": sha256(os.path.join(os.path.dirname(HUMAN_PT_TABLE), "provenance.json")),
            "barcodes_match_object": True, "barcodes_match_v9_object": V9_HUMAN,
            "spearman_with_pt_of_record": float(pd.Series(pt).rank().corr(pd.Series(pt_old).rank())),
            "median_pt_by_label": {str(k): float(v) for k, v in side.groupby("label", observed=True)[HUMAN_PT_COL].median().items()},
            "median_pt_by_timepoint": {str(k): float(v) for k, v in side.groupby("timepoint", observed=True)[HUMAN_PT_COL].median().items()},
        }
    ok = ~np.isnan(pt)
    chk = {"n_obs": int(len(idx)), "n_with_pseudotime": int(ok.sum()), "expected_n_with_pseudotime": EXPECTED_N[sp],
           "pt_min": float(np.nanmin(pt)), "pt_max": float(np.nanmax(pt))}
    if "subtype" in cols:
        c, v = cols["subtype"]
        chk["subtype_column"] = c
        chk["median_pt"] = {s: float(np.median(pt[ok & (v == s)])) for s in ("decidual_boundary", "decidual_proliferative", "decidual_mature") if (ok & (v == s)).any()}
        chk["n_nan_by_subtype"] = {str(s): int(n) for s, n in pd.Series(v[~ok]).value_counts().items()}
    if "day" in cols:
        c, v = cols["day"]
        chk["day_column"] = c
        chk["spearman_pt_day"] = float(pd.Series(pt[ok]).rank().corr(pd.Series(v[ok].astype(float)).rank()))
    auc_path = f"{LANE}/scenic_seeds/{sp}/{AUC}"
    auc = pd.read_csv(auc_path, index_col=0)
    keep = idx[ok]
    miss = keep.difference(auc.index)
    chk["n_with_pseudotime_missing_activity"] = int(len(miss))
    chk["n_activity_nuclei_without_pseudotime"] = int(len(auc.index.difference(keep)))
    assert len(miss) == 0, f"{sp}: {len(miss)} nuclei with pseudotime have no activity score, e.g. {list(miss[:3])}"
    ptk = pd.Series(pt[ok], index=keep)
    A = auc.loc[keep]
    r = spearman(A.to_numpy(), ptk.to_numpy())
    t = pd.DataFrame({"regulon": auc.columns, "rho_pseudotime": r})
    t.to_csv(f"regulon_rho_pseudotime_{sp}{SUFFIX}.tsv", sep="\t", index=False, float_format="%.6g")
    rho[sp] = dict(zip(t.regulon.str[:-3], t.rho_pseudotime))
    chk["n_nuclei_correlated"] = int(len(keep)); chk["n_regulons"] = int(auc.shape[1])
    prov["species"][sp] = chk
    prov["inputs"][f"activity_{sp}"] = {"path": auc_path, "sha256": sha256(auc_path)}
    prov["inputs"][f"object_{sp}"] = {"path": OBJ[sp], "sha256": sha256(OBJ[sp])}
    print(sp, json.dumps(chk), flush=True)

O = pd.read_csv(ORTHO)[["human", "mouse"]].drop_duplicates()
assert O.human.is_unique and O.mouse.is_unique
prov["inputs"]["orthologs"] = {"path": os.path.abspath(ORTHO), "sha256": sha256(ORTHO)}
d = O[O.human.isin(rho["human"]) & O.mouse.isin(rho["mouse"])].rename(columns={"human": "human_tf", "mouse": "mouse_tf"})
d = d.sort_values("human_tf").reset_index(drop=True)
d["rho_human"] = d.human_tf.map(rho["human"])
d["rho_mouse"] = d.mouse_tf.map(rho["mouse"])
d["label"] = d.human_tf.isin(LABELLED).astype(int)
missing = sorted(set(LABELLED) - set(d.human_tf))
d.to_csv(f"regulon_quadrants_input{SUFFIX}.csv", index=False, float_format="%.6g")
prov["n_tf_pairs"] = int(len(d)); prov["labelled"] = sorted(d.human_tf[d.label == 1]); prov["labelled_not_paired"] = missing
json.dump(prov, open(f"scenic_quadrant_rho{SUFFIX}.provenance.json", "w"), indent=1)
print(f"{len(d)} TF pairs with a regulon in both species; {int(d.label.sum())} labelled; not paired: {missing}")
