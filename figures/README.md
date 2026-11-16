# figures: one notebook per figure

| notebook | draws |
|---|---|
| `Figures_Panels` (Çağrı Çevrim) | the cell type panels: UMAPs, dot plots, stacked bars and spatial plots |
| `Fig4_mouse_decidual_maturation` | Fig. 4D: mouse pseudotime on the tissue sections |
| `MI_Genes_Heatmap` (Çağrı Çevrim) | Fig. 5F, from Data S3 (the MI table) and the deposited objects |
| `Fig6_crossspecies_radial_patterning` | Fig. 6: mouse cell type and pseudotime transferred to human by SAMap |
| `S21_S29_pseudotime_trajectories` | Fig. S21 and S29: mouse and human pseudotime trajectories |
| `S23B_notum_vs_pseudotime` | Fig. S23B: *Notum* against pseudotime |
| `S30_mi_gene_selection` | Fig. S30: selection of pseudotime-associated genes by mutual information |
| `S32_samap_supp` | Fig. S32: the SAMap supplement; runs the panel notebooks in `samap/notebooks/5_figures/` |

Fig. 4H and the heatmap of pseudotime-associated genes (Fig. S31) are drawn by `pt_heatmaps/`.

Çağrı Çevrim's notebooks run in his own conda environment; `MI_Genes_Heatmap` takes its file paths in its second
code cell. The others run with the `rsc_uv` kernel: each reads the final objects and writes its panels to a new
folder, both set in its first code cell (`OUT` must be a new folder).
