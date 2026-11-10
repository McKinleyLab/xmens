# figures

One notebook per figure. Each draws its figure's panels, from the final objects and the outputs of the analysis
folders next to this one, into `OUT`:

    PM=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/papermill
    mkdir -p $OUT
    $PM S24_pregnancy_zone_modules.ipynb $OUT/S24_pregnancy_zone_modules.executed.ipynb -p OUT $OUT

The files of record are in `xmens_figures/final_figures/` (synced with Google Drive), one folder per figure, and each
notebook is named after its folder. The notebooks never write there (the setup cell refuses an `OUT` inside
`xmens_figures`): run one, check its panels against the files of record in its last section, then copy them over by
hand.

**Outlines, 8 Nov 2026.** No notebook draws anything yet. Each one lists its panels, the file of record, the code that
drew it and how closely that code's output matches the file of record; the code moves in one panel at a time.

## The notebooks

"Match" compares the file of record with the newest output of the code that drew it: *identical*, byte for byte;
*same drawing*, identical once the SVG date stamp and generated element ids are set aside; *same content*, the same
colours, text and embedded-image pixels (positions not compared); *close*, differences named in the notebook.

| notebook | panels | drawn now by | match |
|---|---|---|---|
| Fig3_mouse_celltypes_spatial | 3A-D, 2 options | xmens_final `30_fig_mouse_tissue.smk` | identical |
| Fig4_mouse_decidual_maturation | 4H | `xmens_pt_heatmaps/` step 3 | identical (unlabelled); same content (labelled, colour bar) |
| Fig5_human_mouse_shared_programs | 5B, 5C, 5E, 5F, 1 option | xmens_final `32_fig_human_tissue.smk` (5B, 5C); `scenic/` step 6 (5E); `xmens_pt_heatmaps/` step 3 (5F) | identical; 5F redrawn with 156 px differing by at most 3 of 255 |
| Fig6_crossspecies_radial_patterning | 6A-E | xmens_final `34_fig_cross_species.smk` (6A, 6D); `samap/` step 5 (6B, 6C, 6E) | identical (6A, 6D); same content (6B, 6E); 6C: same pixels, one legend label changed since |
| S12_mouse_qc_annotation | S12A-D | xmens_final `30_fig_mouse_tissue.smk` | identical |
| S13_mouse_spatial_celltypes | S13 | xmens_final `30_fig_mouse_tissue.smk` | identical |
| S15_mouse_cxcl12 | S15A, S15B | xmens_final `31_fig_mouse_stroma.smk` | same content (A); identical (B) |
| S16_mouse_spatial_fibdec_subtypes | S16 | xmens_final `31_fig_mouse_stroma.smk` | identical |
| S17_mouse_fibdec_subtypes | S17A-C, 1 option | xmens_final `31_fig_mouse_stroma.smk` | same content or same drawing; option identical |
| S20_mouse_pseudotime | composite | `tmp/2026_11_02_pt_supplements/make_pseudotime_supplements_v2.py` | same content |
| S22_mouse_notum | sections; Notum against pseudotime (2 files) | `tmp/2026_11_02_s22_notum_recreate/`; values only in `tmp/2026_11_03_notum_current_objects/` | close (sections: one key label differs); same drawing (dec1); none found (all nuclei) |
| S24_pregnancy_zone_modules | S24A, S24B | `tmp/zone_final_objects_2026-11-02/drive_bundle/final_object_2026-11-02/reproduce_S24.sh` | same drawing; `zone_modules/` draws a newer style |
| S25_human_qc_annotation | S25A-D | xmens_final `32_fig_human_tissue.smk` | identical |
| S26_human_spatial_celltypes | S26A, 8 options for B-C | xmens_final `32_fig_human_tissue.smk`, `33_fig_human_epithelium.smk` | identical |
| S27_human_pseudotime | composite | `tmp/2026_11_02_pt_supplements/make_pseudotime_supplements_v2.py` | same content |
| S28_scenic_supp | S28 | `scenic/` step 6 | identical |
| S29_mi_gene_selection | S29A-D | `tmp/2026_10_30_s29_5f_bigheat/mi_gene_selection_figS29_fibdec_latent/` notebook | same drawing; B identical |
| S30_samap_supp | composite | `samap/` step 5 | same content |
| big_mi_heatmaps | 1 | `xmens_pt_heatmaps/` step 3 | same colours and text; two embedded images 1 px wider |
| C7 | 1 | not found (v7: `xmens_figures/c7_figure/`) | none found |
| IVD_comparison_bars | 20 files | not found (numbers: xmens_final `39_ivd_comparisons.smk`) | none found |
| glands | 1 | `tmp/2026_11_02_epithelial_coarse_clusters/from_final_object/` | close (same colours and text; up to 1,391 px differ) |
| TableM1_marker_dotplots | dot plots, Leiden review | xmens_final `36_fig_marker_dotplots.smk` | not checked (not in final_figures) |
| spatial_qc_supp | 1 | `xmens_figures/spatial_qc_supp/fig_spatial_qc.executed.ipynb` | not checked (not in final_figures) |

`final_figures/bi_mi_gene_heatmaps/` holds an identical copy of the big MI heatmap. Tables are not covered here:
`final_figures/DESeq2_DE_and_IVD_comparisons/` holds two Excel files.

## Not here yet

- No file in final_figures: Fig. 1 and 2, panels 4A-G, 5A and 5D, and Figs. S1-S11, S14, S18, S19, S21 and S23.
  Any of these drawn by code still needs a notebook.
- Planned on the to-do list, no file yet: decidual-subtype lineages from moscot without 0.5 dpi, drawn as origin
  probabilities; the GRN network; functionalis and basalis markers (Gargett lab) with C7; SAMap transfer of the
  decidual subtype labels.

## Conventions

- One notebook per `final_figures` folder, with the folder's name. One SVG per panel in `OUT`, named as its file of
  record.
- Inputs are the final objects and the outputs of the folders in `xmens_scripts`, not `tmp/` or `xmens_final/`.
- The first code cell holds the parameters (`OUT`); the last section compares `OUT` with the files of record.
