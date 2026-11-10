# record/: the pseudotime heatmaps of record (copy of record, 2026-11-07)

Unedited copies from the heatmap package of 2 Nov 2026 (`2026-11-02_pseudotime_heatmaps_final_objects.tar.gz`), which
drew the paper's three pseudotime heatmaps from the final objects. Nothing here is run; the step scripts one folder up
redo the numbers with Monocle 2 itself and the drawing with these renderers cut down to the three figures, and
`compare_with_record` checks their outputs against these files.

| Figure | Figure of record (Drive `xmens_figures/final_figures/`) | Copy here (`package/`) | Made by (`package/code/`) |
|---|---|---|---|
| Fig. 4H | `Fig4_mouse_decidual_maturation/`: `Fig4H_mouse_pseudotime_gene_heatmap_unlabelled_30x59p887mm_not_v9.svg`, `Fig4H_mouse_pseudotime_gene_heatmap_labelled_30x59p887mm_not_v9.svg`, `Fig4H_heatmap_colorbar_30x59p887mm_not_v9.svg` | `yang_induced/`: `heatmap_monocle_induced_z3_box_30x59p887mm.svg`, `heatmap_monocle_induced_labelled_30x59p887mm.svg`, `heatmap_monocle_z3_colorbar_5pt_30x59p887mm.svg`; z matrix `heatmap_monocle_induced_z.tsv` | `yang/monocle_heatmap_repro.R` (curves and z), `yang/render_heatmap_box_svg.py`, `yang/compose_heatmaps_svg.py` |
| Fig. 5F | `Fig5_human_mouse_shared_programs/Fig5F_shared_pseudotime_genes_heatmap_human_mouse_final_objects.svg` | `crossspecies/crossspecies28_monocle_clip3_final_objects_nostrips.svg`; z matrices `..._z_human.csv`, `..._z_mouse.csv` | `crossspecies/monocle_conserved_smooth.R`, `crossspecies/render_crossspecies_heatmaps_nostrips.py` |
| Big MI heatmaps | `big_mi_heatmaps/mi_top_genes_heatmap_centroid_order_final_objects.svg`, gene-order tables in `big_mi_heatmaps/code/` | `mi_heatmaps/mi_top_genes_heatmap_centroid_order.svg`; `..._gene_order_human.csv`, `..._gene_order_mouse.csv` | `crossspecies/monocle_conserved_smooth.R`, `mi/render_full_heatmaps.py` |

Each Drive file is byte-identical to its copy here (sha256 in `SOURCES.tsv`; the Drive folder is mirrored on the
cluster at `/n/lab_storage/eddy_lab/users/nhilgert/xmens_figures/final_figures/`).

How the record was made (`package/README.md`, `package/run_all.sh`):
- `code/inputs/export_heatmap_inputs_final_objects.py` ran on the cluster on 2 Nov 2026 and wrote `inputs/`: per
  nucleus the pseudotime, label, timepoint and total counts, the counts of the genes drawn, and one dispersion table
  per species (`inputs/checks.json` holds the input sha256 and the checks against the Methods).
- `run_all.sh` ran on a Mac the same day: the R scripts made the curves and the Python scripts drew the figures
  (matplotlib 3.11.1, named in the Fig. 5F and MI SVGs; text measured with Helvetica Neue).
- The R scripts are a port of Monocle 2's smoothing written with VGAM, from the monocle-release GitHub master (commit
  7df1050, 2019); they do not call the monocle package. Size factors are total counts over their geometric mean, over
  all nuclei with a pseudotime; the dispersion trend is Monocle's (a + b / mean, fitted with a Gamma GLM and refitted
  without Cook's-distance outliers); each gene is fitted as in `genSmoothCurves()` (VGAM `vglm` of the counts over the
  size factor, rounded, on `sm.ns(Pseudotime, df = 3)`, negative binomial with the dispersion from the trend, epsilon
  0.1) at 100 points; curves become z as in `plot_pseudotime_heatmap()` (log10(curve + 1), z per gene with the n - 1
  SD, clipped at ±3).
- The dispersion tables (`inputs/common/disp_table_*.tsv.gz`) are Monocle's moment estimate with the sample variance
  (n - 1); Monocle uses the population variance (n). `port_dispersion_trend.tsv` holds the trend coefficients the R
  scripts printed in that run, to 4 significant digits; the logs were not kept, and the lines are copied from the
  analysis session's transcript.

Not copied: the per-nucleus counts of the MI-heatmap genes (`inputs/mi/percell_{human,mouse}_topdecile.tsv.gz`, 36 MB;
step 1 writes every gene's counts), the peak-ordered MI heatmaps (not in the paper), and the R scripts' curves
(`work/`, not in the package; their sha256 are in the `*_provenance.json` files).

`SOURCES.tsv`: every file here, its source, its sha256, and the Drive file with the same bytes where there is one.
