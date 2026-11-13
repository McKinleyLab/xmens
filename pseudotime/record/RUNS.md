## 2026-11-12 — human pseudotime, final: root rule restricted to eFib1 (full run from the all-cells object)
Purpose:  final human pseudotime: the 7 Nov pipeline with the root rule restricted to eFib1 (12 Nov meeting decision)
Code:     pseudotime/pseudotime_human_final_1112.ipynb @ 0b83615
Run:      JupyterLab, python /n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python, node holygpu7c0530.rc.fas.harvard.edu, GPU NVIDIA A40,
          started 2026-11-12 19:46, output written 2026-11-12T20:07:27
Settings: ROOT_RULE 'record', ROOT_LABEL 'eFib1', OT_ON_CPU True,
          FIT_LATENT_ON_ALL_NUCLEI False, SCVI_MAX_EPOCHS None,
          SCVI_MODEL_DIR None (latent saved in the output as obsm['X_scVI_pseudotime'])
Inputs:   INPUT_OBJECT     /n/lab_storage/eddy_lab/Lab/mckinley/data/20261107_Human_allcells_final.h5ad  sha256 8d4f0764866d089789ee57e3c3aa6c27a9cd319b06586b9fe4e88d6d09b09b67
          RECORD_LATENT    /n/lab_storage/eddy_lab/users/nhilgert/tmp/2026_10_30_human_fibdec_latent/run4_a40/fibdec/latent_human_stromal_fibdec.h5ad  (HVG check only)
          REFERENCE_OBJECT /n/netscratch/eddy_lab/Everyone/nhilgert/final_objects/Cevrim_Human_Only_Fib_Dec_Cells.h5ad  (deposited object; end comparison only)
Outputs:  /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_human_reroot/20261112_Human_fibdec_final_nick_pdt_full.h5ad  sha256 c7d4bacb40f6472000c372cd005c6712a9c3f0e112ae5b23fb7df04b31dea82b
          executed notebook: ____
Printed:  HVGs identical to record False | scVI epochs 178 |
          moscot converged True on cpu | trajectory DC 1 |
          root GGCACGTTCAGTCCGG-1-secretory_mck6a (eFib1; rule: record, restricted to eFib1) | rho with timepoint 0.269
Embedding vs 7 Nov: ____   (identical diffusion map = the root is the only change)
Status:   final candidate

## 2026-11-12 — human pseudotime, rerun of the final (reproducibility check)
Code:     pseudotime/pseudotime_human_final_1112_rerun.ipynb, a copy of pseudotime/pseudotime_human_final_1112.ipynb (18 vs 18 code cells; differing code cells [0])
Output:   /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_human_reroot_rerun/20261112_Human_fibdec_final_nick_pdt_full.h5ad  sha256 bc16b54469b15b2915e41a7bfe524ea4cf1d00baf97daacce36a7e724f04a82c  (written 2026-11-12T21:12:38)
Check:    h5diff vs /n/lab_storage/eddy_lab/Lab/mckinley/data/20261112_Human_fibdec_final_nick_pdt.h5ad, excluding /uns/pseudotime/computed: identical (h5diff exit 0)
Executed: pseudotime/record/pseudotime_human_final_1112_rerun.executed_20261112.ipynb
Status:   check

## 2026-11-12 — mouse pseudotime, rerun of the 10 Nov object of record (reproducibility check)
Code:     pseudotime/pseudotime_mouse_final_1110_rerun.ipynb, a copy of pseudotime/pseudotime_mouse_final_1110.ipynb (18 vs 18 code cells; differing code cells [0])
Output:   /n/lab_storage/eddy_lab/Lab/mckinley/data/20261110_mouse_only_eFib_dec_final_nick_pdt_rerun.h5ad  sha256 126024141c6dcb96145c8a737cb2fcdd31c648dde322531824cd1600f8fb7dbf  (written 2026-11-12T21:18:30)
Check:    h5diff vs /n/lab_storage/eddy_lab/Lab/mckinley/data/20261110_mouse_only_eFib_dec_final_nick_pdt.h5ad (sha256 as pinned), excluding /uns/pseudotime/computed: identical (h5diff exit 0)
Executed: pseudotime/record/pseudotime_mouse_final_1110_rerun.executed_20261112.ipynb
Status:   check
