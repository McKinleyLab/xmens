#!/usr/bin/env python
"""Fig. 5E, step 1: Spearman rho of every floor-6 consensus regulon's activity with the pseudotime, both species, and the
quadrant-plot input table (rule programs_regulon_table, 35_fig_programs.smk).

Port of scenic_quadrant_rho_fibdec_latent.py (2026-10-30, sha256 e2ea80fe...), the script that made the table of record
(resources/figures/panels_programs/regulon_quadrants_input_fibdec_latent.csv): same computation, every path an argument.
  activity   the CLI AUCell matrix of the floor-6 consensus (auc_mtx_weighted.csv) per species
  pseudotime obs[--pt-col] of each species' final fibroblast/decidual object, the rank-normalized pseudotime of record
             (mouse: the 22 Oct pt_harmonized rank, NaN for the apoptotic decidual nuclei, which are dropped; human: the
             30 Oct refit on the fibroblast/decidual latent); nothing is re-ranked
  orthologs  Ensembl one-to-one human-mouse symbol pairs (rule mi_overlap)
Spearman rho = Pearson correlation of average ranks; activity cast to float32 before ranking. Nuclei are matched to the
activity matrix by obs_names; the script stops unless every nucleus with a pseudotime has an activity score.
Outputs: regulon_rho_pseudotime_{mouse,human}.tsv, regulon_quadrants_input.csv (human_tf, mouse_tf, rho_human, rho_mouse,
label), provenance.json (sha256 of every input and the sanity checks).
"""
import argparse, datetime, hashlib, json, os
import h5py, numpy as np, pandas as pd
from scipy.stats import rankdata

# Subtype labels serve only the sanity checks in the provenance (median pseudotime and nuclei without pseudotime per
# subtype). obs is read by column name (index, pseudotime, subtype) and no other obs column is read.
SUBTYPE_COL = "cell_subtype"


def sha256(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
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


def read_obs(path, pt_col):
    """obs_names, the pseudotime column and the subtype column."""
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        idx = strs(obs[obs.attrs["_index"]])
        pt = obs_col(obs, pt_col).astype(float)
        subtype = obs_col(obs, SUBTYPE_COL)
    return pd.Index(idx), pt, subtype


def spearman(A, pt):
    """Spearman rho of each column of A (nuclei x regulons) with pt."""
    pr = rankdata(pt); pr -= pr.mean()
    RA = rankdata(np.asarray(A, dtype=np.float32), axis=0); RA -= RA.mean(axis=0)
    return (pr @ RA) / (np.sqrt((RA * RA).sum(axis=0)) * np.sqrt((pr * pr).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--activity-mouse", required=True); ap.add_argument("--activity-human", required=True)
    ap.add_argument("--object-mouse", required=True); ap.add_argument("--object-human", required=True)
    ap.add_argument("--pt-col", default="pseudotime", help="obs column of both objects: the pseudotime of record")
    ap.add_argument("--orthologs", required=True, help="ortholog_pairs_one2one.csv with columns human, mouse")
    ap.add_argument("--labelled", required=True, help="comma-separated human TF symbols flagged in the label column")
    ap.add_argument("--expected-n", default="", help="e.g. mouse=42290,human=47533 (recorded, not enforced)")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    labelled = [x for x in a.labelled.split(",") if x]
    expected = {k: int(v) for k, v in (kv.split("=") for kv in a.expected_n.split(",") if kv)}
    activity = {"mouse": a.activity_mouse, "human": a.activity_human}
    objects = {"mouse": a.object_mouse, "human": a.object_human}

    prov = {"written": datetime.datetime.now().isoformat(timespec="seconds"),
            "pseudotime_column": {"mouse": a.pt_col, "human": a.pt_col},
            "inputs": {}, "species": {}}
    rho = {}
    for sp in ("mouse", "human"):
        idx, pt, subtype = read_obs(objects[sp], a.pt_col)
        ok = ~np.isnan(pt)
        chk = {"n_obs": int(len(idx)), "n_with_pseudotime": int(ok.sum()), "expected_n_with_pseudotime": expected.get(sp),
               "pt_min": float(np.nanmin(pt)), "pt_max": float(np.nanmax(pt))}

        # Decidual subtypes should sit late in pseudotime, and in mouse only the apoptotic subtype lacks a value.
        chk["subtype_column"] = SUBTYPE_COL
        median_by_subtype = pd.Series(pt[ok]).groupby(subtype[ok]).median()
        chk["median_pt"] = {str(s): float(v) for s, v in median_by_subtype.items()}
        chk["n_nan_by_subtype"] = {str(s): int(n) for s, n in pd.Series(subtype[~ok]).value_counts().items()}

        auc = pd.read_csv(activity[sp], index_col=0)
        keep = idx[ok]
        miss = keep.difference(auc.index)
        chk["n_with_pseudotime_missing_activity"] = int(len(miss))
        chk["n_activity_nuclei_without_pseudotime"] = int(len(auc.index.difference(keep)))
        if len(miss):
            raise SystemExit(f"{sp}: {len(miss)} nuclei with pseudotime have no activity score, e.g. {list(miss[:3])}")
        ptk = pd.Series(pt[ok], index=keep)
        A = auc.loc[keep]
        r = spearman(A.to_numpy(), ptk.to_numpy())
        t = pd.DataFrame({"regulon": auc.columns, "rho_pseudotime": r})
        t.to_csv(os.path.join(a.out_dir, f"regulon_rho_pseudotime_{sp}.tsv"), sep="\t", index=False, float_format="%.6g")
        rho[sp] = dict(zip(t.regulon.str[:-3], t.rho_pseudotime))
        chk["n_nuclei_correlated"] = int(len(keep)); chk["n_regulons"] = int(auc.shape[1])
        prov["species"][sp] = chk
        prov["inputs"][f"activity_{sp}"] = {"path": activity[sp], "sha256": sha256(activity[sp])}
        prov["inputs"][f"object_{sp}"] = {"path": objects[sp], "sha256": sha256(objects[sp])}
        print(sp, json.dumps(chk), flush=True)

    O = pd.read_csv(a.orthologs)[["human", "mouse"]].drop_duplicates()
    if not (O.human.is_unique and O.mouse.is_unique):
        raise SystemExit(f"{a.orthologs}: not one-to-one")
    prov["inputs"]["orthologs"] = {"path": os.path.abspath(a.orthologs), "sha256": sha256(a.orthologs)}
    d = O[O.human.isin(rho["human"]) & O.mouse.isin(rho["mouse"])].rename(columns={"human": "human_tf", "mouse": "mouse_tf"})
    d = d.sort_values("human_tf").reset_index(drop=True)
    d["rho_human"] = d.human_tf.map(rho["human"])
    d["rho_mouse"] = d.mouse_tf.map(rho["mouse"])
    d["label"] = d.human_tf.isin(labelled).astype(int)
    missing = sorted(set(labelled) - set(d.human_tf))
    d.to_csv(os.path.join(a.out_dir, "regulon_quadrants_input.csv"), index=False, float_format="%.6g")
    prov["n_tf_pairs"] = int(len(d)); prov["labelled"] = sorted(d.human_tf[d.label == 1]); prov["labelled_not_paired"] = missing
    json.dump(prov, open(os.path.join(a.out_dir, "provenance.json"), "w"), indent=1)
    print(f"{len(d)} TF pairs with a regulon in both species; {int(d.label.sum())} labelled; not paired: {missing}", flush=True)


if __name__ == "__main__":
    main()
