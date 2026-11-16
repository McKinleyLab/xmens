# notebooks

## The deposited fibroblast/decidual objects (14 Nov 2026)

Run in this order, in JupyterLab with the `rsc_uv` kernel. Each notebook only reads its inputs and stops rather than
overwrite an output.

1. `stro7_score_transfer_20261114.ipynb`: the Stro-7 module score (`Ulrich2024_sg_Stro7_top50`), transferred from the
   human to the mouse nuclei through the SAMap graph; writes `stro7_top50_20261114/`.
2. `add_columns_20261114.ipynb`: adds pseudotime (`pseudotime/`), the SAMap transfers (`samap/`) and the Stro-7
   score to the annotated objects of 11 Nov.
3. `add_pt_samap_pt_info_fibdec_20261114.ipynb`: writes the SAMap object with labels and pseudotime from the
   fibroblast/decidual objects (`20261114_samap_final_ensembl_ids_annotated.h5ad`, `..._pseudotime.h5ad`).
4. `update_cevrim_objects_20261114.ipynb`: replaces the three SAMap-derived columns of the deposited objects with
   those of the 14 Nov SAMap run and writes the deposited objects.

## Other analyses

- `functionalis_basalis_markers/`: genes higher in decidual than in fibroblast nuclei, or the reverse, in both
  species (pseudobulk DESeq2), plotted on every section; run by its `.sbatch`.
- Çağrı Çevrim's notebooks: `Integration, Cleaning, and Annotation.ipynb` (step 2), `Figures_Panels.ipynb`,
  `Human-Stro7_scoring.ipynb`, `Mouse DEGs.ipynb`, `CC3_pH3_Segmentation.ipynb` and the bulk RNA-seq DEG analyses
  `01`-`04_DEG_Analysis_*.ipynb` (environment: `environments/BulkRnaSeq.yml`).
