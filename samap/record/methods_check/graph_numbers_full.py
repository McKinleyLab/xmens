"""Independent check of numbers the Methods report from the SAMap joint graph (written while verifying the Methods,
7 Nov 2026; not part of the pipeline): SAMap's alignment scores between compartments, recomputed from the graph with
the formula of samap.analysis.scores._compute_csim (n_top=0), and per compartment the total cross-species edge weight
and number of cross-species edges per nucleus. Validated on two graphs for which SAMap's own get_mapping_scores output
exists: the fit of record with its fit-time labels (record mapping_celltype.csv) and the fresh fit of
2026_11_06_xmens_samap with the final labels (its 2_samap/mapping_celltype.csv). Then applied to the fit of record with
the final labels (Samap_Final.h5ad), which gives the values in the Methods."""

import json, sys
import numpy as np, pandas as pd, scipy.sparse as sp, h5py
from anndata.io import read_elem

FO = "/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects"
REC = "/n/lab_storage/eddy_lab/users/nhilgert/code/samap-scvi/full_fibdec/runs_raw/vanilla_raw"
RUNF = "/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_xmens_samap"
out = {}

def load(p, cols):
    with h5py.File(p, "r") as f:
        obs = pd.DataFrame({c: read_elem(f["obs"][c]) for c in cols})
        G = read_elem(f["obsp/connectivities"]).tocsr()
        K = float(read_elem(f["uns/mapping_K"]))
        obs.index = read_elem(f["obs"][f["obs"].attrs["_index"]])
    G.eliminate_zeros()
    return obs, G, K

def analyse(name, obs, G, K, lab):
    spc = obs["species"].astype(str).to_numpy()
    lab = np.asarray(lab).astype(str)
    cl = np.char.add(np.char.add(spc, "_"), lab)
    ix = np.unique(spc, return_index=True)[1]; skeys = spc[np.sort(ix)]
    clu = np.concatenate([np.unique(cl[spc == s]) for s in skeys])
    coo = G.tocoo(); f = spc[coo.row] != spc[coo.col]
    r_, c_, w = coo.row[f], coo.col[f], coo.data[f]
    j = pd.Index(clu).get_indexer(cl[c_])
    M = sp.coo_matrix((w, (r_, j)), shape=(G.shape[0], clu.size)).toarray()   # cell x cluster cross-species weight
    CS = np.vstack([M[cl == c].mean(0) for c in clu])                        # directional means
    A = np.maximum(CS, CS.T) / K                                               # SAMap alignment score
    res = {"K": K, "symmetric_graph": bool(abs(G - G.T).sum() < 1e-6), "clusters": clu.tolist(),
           "alignment": {f"{a}|{b}": float(A[i, k]) for i, a in enumerate(clu) for k, b in enumerate(clu) if a.split("_")[0] != b.split("_")[0] and i < k},
           "directional_mean_over_K": {f"{a}->{b}": float(CS[i, k] / K) for i, a in enumerate(clu) for k, b in enumerate(clu) if a.split("_")[0] != b.split("_")[0]}}
    tot = M.sum(1); nnz = np.bincount(r_, minlength=G.shape[0])
    comp = {}
    for c in clu:
        m = cl == c
        comp[c] = {"n": int(m.sum()), "mean_total_weight": float(tot[m].mean()),
                   "max_total_weight": float(tot[m].max()), "n_total_gt_K": int((tot[m] > K + 1e-9).sum()),
                   "cross_neighbors_min_median_max": [int(nnz[m].min()), float(np.median(nnz[m])), int(nnz[m].max())],
                   "n_no_cross_neighbor": int((tot[m] == 0).sum())}
    res["compartments"] = comp
    for s in skeys:
        m = spc == s
        res[f"{s}_no_cross_neighbor"] = int((tot[m] == 0).sum())
        res[f"{s}_n"] = int(m.sum())
    out[name] = res
    print(name, json.dumps(res["alignment"]), flush=True)

# 1) record fit, labels at fit time (validates the formula against the record's SAMap mapping_celltype.csv)
obs, G, K = load(f"{REC}/samap_transferred.h5ad", ["species", "cell_type_annotations", "hs_cell_type_annotations", "mm_cell_type_annotations"])
print("record joint: cell_type_annotations nulls", int(obs["cell_type_annotations"].isna().sum()), flush=True)
lab = obs["cell_type_annotations"].astype(object)
analyse("record_fit_fit_time_labels", obs, G, K, lab)
out["record_mapping_celltype_csv"] = pd.read_csv(f"{REC}/mapping_celltype.csv", index_col=0).to_dict()
# 2) record fit, final labels (Samap_Final of record)
obs2, G2, K2 = load(f"{FO}/Samap_Final.h5ad", ["species", "cell_type"])
assert (obs2.index == obs.index).all()
out["graph_identical_record_joint_vs_Samap_Final"] = bool((G != G2).nnz == 0)
analyse("record_fit_final_labels", obs2, G2, K2, obs2["cell_type"].astype(object))
# 3) fresh fit, final labels (validates against SAMap's own get_mapping_scores in that run)
obs3, G3, K3 = load(f"{RUNF}/3_transfer/Samap_Final.h5ad", ["species", "cell_type"])
analyse("fresh_fit_final_labels", obs3, G3, K3, obs3["cell_type"].astype(object))
out["fresh_mapping_celltype_csv"] = pd.read_csv(f"{RUNF}/2_samap/mapping_celltype.csv", index_col=0).to_dict()
json.dump(out, open(sys.argv[1], "w"), indent=1)
print("done")

