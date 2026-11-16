# figures: one notebook per figure

Each notebook reads the final objects and writes its panels to a new folder; both are set in its first code cell
(`OUT` must not exist yet, so nothing is overwritten). Run all cells with the `rsc_uv` kernel.

| notebook | draws |
|---|---|
| `Fig4_mouse_decidual_maturation` | Fig. 4D: mouse pseudotime on the tissue sections |
| `Fig6_crossspecies_radial_patterning` | Fig. 6: mouse cell type and pseudotime transferred to human by SAMap |
| `S21_S29_pseudotime_trajectories` | mouse and human pseudotime trajectories |
| `S23B_notum_vs_pseudotime` | *Notum* against pseudotime |
| `S30_mi_gene_selection` | selection of pseudotime-associated genes by mutual information |
| `S32_samap_supp` | Fig. S32: the SAMap supplement; runs the panel notebooks in `samap/notebooks/5_figures/` |

File names carry the figure numbers of an earlier draft; the table gives the current number where it differs.
Fig. 4H and the heatmap of pseudotime-associated genes (Fig. S31) are drawn by `pt_heatmaps/`.

`MI_Genes_Heatmap.ipynb` draws Fig. 5F from the MI table (Data S3) and the deposited fibroblast/decidual objects,
whose paths are set in its second code cell.

`Figures_Panels.ipynb` (Ç. Çevrim) draws the cell type panels: UMAPs, dot plots, stacked bars and spatial plots. It
runs in his own conda environment.
