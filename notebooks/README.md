# notebooks

The deposited fibroblast/decidual objects, made in this order (JupyterLab, `rsc_uv` kernel):

1. `stro7_score_transfer_20261114.ipynb`: the Stro-7 score, transferred from human to mouse through the SAMap graph
2. `add_columns_20261114.ipynb`: adds pseudotime, the SAMap transfers and the Stro-7 score
3. `add_pt_samap_pt_info_fibdec_20261114.ipynb`: the SAMap object with labels and pseudotime
4. `update_cevrim_objects_20261114.ipynb`: the columns from the 14 Nov SAMap run; writes the deposited objects

By Çağrı Çevrim:

- `Integration, Cleaning, and Annotation.ipynb`: integration, cleaning and annotation of the mouse and human data
- `Human-Stro7_scoring.ipynb`: the Stro-7 module score (Ulrich 2024, top 50 genes) of the human nuclei, by subtype
  and on the sections
- `Mouse DEGs.ipynb`: X-Mens decidual vs fibroblast DEGs compared with oil-induced artificial decidualization
- `CC3_pH3_Segmentation.ipynb`: CC3+ and pH3+ cells across the decidualization stage series

They run in his conda environments: `BulkRnaSeq` (`../environments/BulkRnaSeq.yml`) for `Mouse DEGs`, `single-cell`
for `Human-Stro7_scoring`, `pt2.3.0_cuda12.1` for the annotation and `cellpose` for the segmentation.
