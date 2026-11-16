# pseudotime: diffusion pseudotime of the fibroblast and decidual nuclei

One notebook per species, run in JupyterLab with the `rsc_uv` kernel on a GPU node (scVI). Each reads an annotated
object and writes the same object with `obs['pseudotime']`. Paths and settings are in the first code cell, and the
notebook stops rather than overwrite its output.

| notebook | reads | object of record |
|---|---|---|
| `pseudotime_mouse_final_1110.ipynb` | `20261110_mouse_only_eFib_dec_final.h5ad` | `20261110_mouse_only_eFib_dec_final_nick_pdt.h5ad` |
| `pseudotime_human_final_1112.ipynb` | `20261107_Human_allcells_final.h5ad` | `20261112_Human_fibdec_final_nick_pdt.h5ad` |

Method, as in the Methods ("Pseudotime"): a stromal scVI latent space (3,000 HVGs, 30 dimensions, seed 0); moscot
couplings between timepoints (mouse) or from secretory to menstrual (human); a CellRank RealTimeKernel; a
15-component diffusion map of its transition matrix; diffusion pseudotime from a root fibroblast, rank-normalized.
Apoptotic decidual nuclei (decAp) get no pseudotime. The human root is restricted to the outer fibroblasts (eFib1).

The human root: `pseudotime_human_final_1107.ipynb` is the human run of 7 Nov 2026, whose root rule picked an eFib2
nucleus. `pseudotime_human_reroot_1112.ipynb` applies the eFib1 root to that run and changes nothing else (it reads
the 7 Nov object); `pseudotime_human_final_1112.ipynb` then runs the whole pipeline with the eFib1 root. Only the last
made an object used in the paper.

The notebooks are saved as run, with their outputs. The runs, their reruns (identical objects) and the runs' logs are
in commit `3605739`.
