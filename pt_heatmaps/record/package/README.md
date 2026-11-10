# Pseudotime heatmaps from the final objects (2 Nov 2026)

All pseudotime heatmaps, regenerated from the final objects in
`/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects/` (`obs['pseudotime']`) and the MI of record, with the code
and inputs to rebuild them.

| folder | figure |
|---|---|
| `yang_induced/` | Yang et al. (2023) 28-gene panel, mouse nuclei from induced timepoints (n = 32,957), Monocle 2 smoothing, z clipped at ±3. Box (30 × 59.887 mm), colorbars (3.6 pt and 5 pt labels), labelled panel and PNG preview. |
| `crossspecies/` | 28 conserved one-to-one orthologs, human (left) and mouse (right), same row order, no composition strips. SVG/PDF/PNG, z matrices, gene list. |
| `mi_heatmaps/` | Every top pseudotime-associated gene in the one-to-one-ortholog background (1,201 human, 1,133 mouse; 373 in both species' top sets, flagged in the `shared` column of the gene-order tables; `code/mi/render_full_heatmaps.py --shared-column` adds the black/white marker column, omitted here), rows ordered by peak pseudotime (`peak_order`) or by pseudotime centroid (`centroid_order`). SVG/PDF/PNG and gene-order tables. |

## Inputs

- Objects: `/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad` (sha256 573c1f8a6a384ccc…) and `/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects/Cevrim_Human_Only_Fib_Dec_Cells.h5ad` (sha256 00dcc9686de56b43…). `obs['pseudotime']` is the
  rank-normalized pseudotime of record. Mouse: identical to `pt_rank` of `/n/netscratch/eddy_lab/Everyone/nhilgert/2026_10_22_harmonized_pt/runs/mouse_pca_released/cells.parquet` (sha256 f1160f48c8a052a8…)
  (max |difference| 0.0, same 42,290 nuclei; the 907 decAp nuclei carry none).
  Human: identical to the 30 Oct 2026 refit, `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/interim/human_pt_2026_10_30_run4_a40/cells.parquet` (sha256 956a07cdb97d9a99…) (max |difference| 0.0, 47,533 nuclei).
- MI of record: `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/results/mi_dpt/pt_harmonized_mouse_nall/mi_all.csv` (sha256 7351f3fe7729e756…) and `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/results/mi_dpt/pt_harmonized_human_nall/mi_all.csv` (sha256 6bcfc0a3fc47fc9c…); ortholog pairs `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/results/mi_dpt/overlap_of_record/ortholog_pairs_one2one.csv` (sha256 8e00db6a2346d92f…).
- `inputs/` was written on the cluster by `code/inputs/export_heatmap_inputs_final_objects.py` (Slurm job
  e365ddba, 2 Nov 2026); `inputs/checks.json` holds every check below and the sha256 of every input.

## Agreement with the Methods (checked on the final objects)

| Methods statement | final objects |
|---|---|
| pseudotime vs timepoint, ρ = 0.62 (mouse), 0.32 (human), Methods day encoding | 0.6219, 0.3198 |
| MI: 14,261 / 15,008 genes tested; 14,261 / 14,887 significant; top 1,426 / 1,489 | 14,261 / 15,008; 14,261 / 14,887; 1,426 / 1,489 |
| top-set thresholds MI > 0.102 / 0.039 bits | 0.1015 / 0.0394 bits |
| background N = 9,587; K = 1,133; n = 1,201; k = 373; expected 141.9; 2.63-fold; p = 2.0 × 10⁻⁸⁴ | 9,587; 1,133; 1,201; 373; 141.9; 2.63-fold; p = 2.0 × 10⁻⁸⁴ (log10 p = -83.7099) (= overlap_of_record: True) |
| cross-species pool: 273 of the 9,638 orthologs tested in both, a subset of the 373 | 273 of 9,638; 273 of 273 in the 373 |
| 92 genes with r > 0.8; 14 early + 14 late shown, r = 0.88–0.99, ordered by centroid | 92 ({'late': 62, 'early': 29, 'mid': 1}); 28 shown, r = 0.882–0.991; same genes and order as the pipeline's `overlap_of_record` selection: True |
| Monocle dispersion trend over the 27,064 genes detected in > 1 nucleus | mouse 27,064 (human, cross-species and MI heatmaps: 31,735) |
| smoothing on all fibroblast/decidual nuclei, 42,290 mouse / 47,533 human | 42,290 / 47,533 |
| Yang panel: induced timepoints, re-ranked within the subset | n = 32,957 |

## Changes from the previous Drive versions

- `yang_induced/` is byte-identical to `2026-10-29_pseudotime_heatmaps_final/yang_induced/` (mouse pseudotime and
  counts are unchanged; only the timepoint labels in the input table are now the +0.5 dpi names).
- `crossspecies/`: the human panel is on the refit human pseudotime. The selection rule picks two different genes:
  AR (early) and LAMP1 (late) drop out, PRKCE (early) and P4HB (late) come in; RGMA and MATN2 swap rows 1 and 2.
  Rows of the 26 kept genes: mouse unchanged; human per-gene r with the old rows ≥ 0.994.
- `mi_heatmaps/`: human set 1,191 → 1,201 genes (120 out, 130 in), mouse 1,132 → 1,133 (1 in), shared 371 → 373
  (supersedes the centroid-ordered `big_migenes_full_{human,mouse}.svg` in `xmens_figures/big_migenes_heatmaps/`).
- Dispersion tables use Monocle's estimator with the sample variance (n − 1) for both species and all figures (the
  previous cross-species and MI tables used n; effect on the mouse cross-species rows < 0.001 in z).

## Regenerate

    bash run_all.sh [NCORES] [all|yang|crossspecies|mi]

Requires R with VGAM, splines, parallel, pheatmap and colorRamps, and Python 3 with numpy, pandas and matplotlib
(cairosvg optional, for the Yang PNG preview). The whole run takes about a minute with 14 cores.
To rebuild `inputs/` from the objects, run on the cluster (Slurm, ~6 min, 96 GB):

    python code/inputs/export_heatmap_inputs_final_objects.py \
      --mouse-h5ad <final_objects>/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad --human-h5ad <final_objects>/Cevrim_Human_Only_Fib_Dec_Cells.h5ad \
      --mi-dir xmens_final/results/mi_dpt --pairs xmens_final/results/mi_dpt/overlap_of_record/ortholog_pairs_one2one.csv \
      --pt-record-mouse <netscratch>/2026_10_22_harmonized_pt/runs/mouse_pca_released/cells.parquet \
      --pt-record-human xmens_final/data/interim/human_pt_2026_10_30_run4_a40/cells.parquet \
      --record-overlap xmens_final/results/mi_dpt/overlap_of_record/methods_overlap.json \
      --record-selection xmens_final/results/mi_dpt/overlap_of_record/crossspecies_heatmap_selection/conserved_profiles_selected_genes.csv \
      --previous-selection xmens_final/results/mi_dpt/overlap_harmonized_nall/crossspecies_heatmap_selection/conserved_profiles_selected_genes.csv \
      --out inputs_new
