# samap: 2026_11_12_xmens_samap_final

Written 2026-11-13 03:47:22 by record/runs_report.py

OUT: /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_xmens_samap_final

## Steps (sacct)

- 1_inputs 60181191: COMPLETED (exit 0:0), 0:02:13, holy8a25503, 2026-11-12 21:54 to 2026-11-12 21:56
- 2_samap 60192327: COMPLETED (exit 0:0), 0:59:10, holy8a27311, 2026-11-12 22:35 to 2026-11-12 23:34
- 3_transfer 60192344: COMPLETED (exit 0:0), 0:01:53, holy8a25505, 2026-11-12 23:35 to 2026-11-12 23:37
- 4_statistics 60192359: COMPLETED (exit 0:0), 0:01:33, holy8a25505, 2026-11-12 23:38 to 2026-11-12 23:39
- 5_figures 60192360: COMPLETED (exit 0:0), 0:05:43, holy8a25505, 2026-11-12 23:40 to 2026-11-12 23:46

## Outputs

107 files, 15.0 GB; sha256 of each in samap/record/manifests/2026_11_12_xmens_samap_final.sha256 (sha256 ef4410acd828d560...)

## Numbers

- mean_alignment: 0.6731769968221357
- seed: 0
- hnswlib_threads: 1

## Compared with /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_xmens_samap_seeded_A

108 files: 37 identical bytes; 19 same content (12 .ipynb, 7 .pdf; time stamps or paths only); 1 not compared; 50 differ; 1 in one folder only

Files not identical in bytes, by folder:

- .: 1 only in other
- 1_inputs: 2 differs
- 2_samap: 2 differs, 1 not compared
- 3_transfer: 2 differs
- 4_statistics: 1 differs
- 4_statistics/graph_block: 1 differs
- 4_statistics/ptnull: 3 differs
- 4_statistics/tables: 2 differs
- 5_figures/fig6: 10 differs, 7 same content
- 5_figures/sfig_samap: 24 differs
- notebooks: 3 differs, 1 same content
- notebooks/5_figures: 11 same content

| file | result | detail |
|---|---|---|
| 1_inputs/hs_fib_dec_all_raw.h5ad | differs | hdf5: 0 objects in one file only, 1 differ: uns/samap_input |
| 1_inputs/inputs.json | differs | json: hs/source: /n/lab_storage/edd..._fibdec_final_nick_pdt.h5ad vs /n/lab_storage/edd..._fibdec_final_nick_pdt.h5ad; hs/source_sha256: c7d4bacb40f6472000...0e112ae5b23fb7df04b31dea82b vs 12429594a2e71ba398...d1f026e11f7ade8856a10a4c425; written: 2026-11-12T21:56:23 vs 2026-11-12T00:34:30 |
| 2_samap/meta.json | differs | json: finished: 2026-11-12T23:33:58 vs 2026-11-12T01:19:44; host: holy8a27311.rc.fas.harvard.edu vs holygpu8a26504.rc.fas.harvard.edu; inputs/hs/sha256: a456275f5a07ff0956...dfc7f3bab52358df1b6f9a7446c vs 6413dc7255dd325457...38c3e14cca8480c40379a62f549; started: 2026-11-12T22:37:52 vs 2026-11-12T00:37:56 |
| 2_samap/samap.pkl | not compared | pickle; not loaded (its contents are compared through the h5ad outputs) |
| 2_samap/samap_joint.h5ad | differs | hdf5: 0 objects in one file only, 1 differ: uns/provenance |
| 3_transfer/Samap_Final.h5ad | differs | hdf5: 0 objects in one file only, 1 differ: uns/provenance |
| 3_transfer/summary.json | differs | json: samap_final/joint_sha256: b48e787da8cab11609...2c283748f0fab4c9ff8252e86eb vs 792dec74eda21fbc7f...af1436422bbda5f36554ec83bb4 |
| 4_statistics/graph_block/meta.json | differs | json: joint_sha256: 7901c4992351a7be8b...624800982ae861d137363860e80 vs 1b4d90f4f0c981a4f3...e016f6438c8d3911240b419f59e |
| 4_statistics/ptnull/null_replicates.csv | differs | text: 2,001 of 2,002 lines differ, first at line 2 |
| 4_statistics/ptnull/null_summary.json | differs | json: by_design/observed/moran_all/sd: nan vs nan; by_design/observed/moran_menstrual/sd: nan vs nan; by_design/observed/moran_secretory/sd: nan vs nan; by_design/observed/moran_within_decidual/sd: nan vs nan; by_design/observed/moran_within_fibroblast/sd: nan vs nan; by_design/observed/rho_all/max: 0.7623 vs 0.6659; by_design/observed/rho_all/mean: 0.7623 vs 0.6659; by_design/observed/rho_all/sd: nan vs nan; ... 29 more |
| 4_statistics/ptnull/ptnull_statistics_table.csv | differs | text: 3 of 7 lines differ, first at line 2 |
| 4_statistics/statistics.json | differs | json: graph_block/joint_sha256: 7901c4992351a7be8b...624800982ae861d137363860e80 vs 1b4d90f4f0c981a4f3...e016f6438c8d3911240b419f59e; ptnull/human_pt_moran/moran_all: 0.6543 vs 0.6319; ptnull/human_pt_moran/moran_menstrual: 0.7199 vs 0.6429; ptnull/human_pt_moran/moran_secretory: 0.5888 vs 0.6209; ptnull/human_pt_moran/moran_within_decidual: 0.6671 vs 0.6197; ptnull/human_pt_moran/moran_within_fibroblast: 0.4956 vs 0.4513; ptnull/table/rho_all/observed: 0.762 vs 0.666; ptnull/table/rho_all/uniform_max: 0.166 vs 0.191; ... 14 more |
| 4_statistics/tables/human_cells.csv | differs | text: 44,882 of 46,819 lines differ, first at line 2 |
| 4_statistics/tables/human_plot_table.parquet | differs | parquet table |
| 5_figures/fig6/20261027-Human-menstrual-transferred_pseudotime-spatial.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/20261027-Human-menstrual-transferred_pseudotime-spatial.svg | differs | svg (date): 28,377 of 42,562 lines differ, first at line 34 |
| 5_figures/fig6/20261027-Human-secretory-transferred_pseudotime-spatial.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/20261027-Human-secretory-transferred_pseudotime-spatial.svg | differs | svg (date): 36,043 of 54,077 lines differ, first at line 34 |
| 5_figures/fig6/20261027-spatial_key_transferred_pseudotime.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/20261027-spatial_key_transferred_pseudotime.svg | differs | svg (date): 5 of 86 lines differ, first at line 41 |
| 5_figures/fig6/fig6_umap_transferred_cell_type.provenance.json | differs | json: git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; host: holy8a25505.rc.fas.harvard.edu vs holy7c10412.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...6cf3cf39a8678c1c042fd10d'}] vs [{'path': '<RUN>/4...6cf3cf39a8678c1c042fd10d'}]; written: 2026-11-12T23:46:05 vs 2026-11-12T01:29:24 |
| 5_figures/fig6/fig6c_umap_transferred_pseudotime.provenance.json | differs | json: git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; host: holy8a25505.rc.fas.harvard.edu vs holy7c10412.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...d78495249e80edfe2cf8006b'}] vs [{'path': '<RUN>/4...528e8f98b1712f277f6a4863'}]; written: 2026-11-12T23:44:49 vs 2026-11-12T01:28:02 |
| 5_figures/fig6/fig6e_sections_transferred_pseudotime.provenance.json | differs | json: git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; host: holy8a25505.rc.fas.harvard.edu vs holy7c10412.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...d78495249e80edfe2cf8006b'}] vs [{'path': '<RUN>/4...528e8f98b1712f277f6a4863'}]; written: 2026-11-12T23:45:20 vs 2026-11-12T01:28:26 |
| 5_figures/fig6/umap_dark_transferred_cell_type.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/umap_dark_transferred_cell_type.svg | differs | svg (date): 11 of 176 lines differ, first at line 41 |
| 5_figures/fig6/umap_dark_transferred_cell_type_humanonly.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/umap_dark_transferred_cell_type_humanonly.svg | differs | svg (date): 8 of 153 lines differ, first at line 43 |
| 5_figures/fig6/umap_dark_transferred_pseudotime.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/umap_dark_transferred_pseudotime.svg | differs | svg (date): 12 of 181 lines differ, first at line 41 |
| 5_figures/fig6/umap_dark_transferred_pseudotime_humanonly.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/umap_dark_transferred_pseudotime_humanonly.svg | differs | svg (date): 9 of 158 lines differ, first at line 43 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_a.provenance.json | differs | json: git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; host: holy8a25505.rc.fas.harvard.edu vs holy7c10412.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...edeab02db39347443335007a'}] vs [{'path': '<RUN>/4...5e45c506e9324572788c4c1e'}]; written: 2026-11-12T23:41:09 vs 2026-11-12T01:25:12 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_a.svg | differs | svg (date): 5 of 125 lines differ, first at line 43 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_b.provenance.json | differs | json: git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; host: holy8a25505.rc.fas.harvard.edu vs holy7c10412.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...edeab02db39347443335007a'}] vs [{'path': '<RUN>/4...5e45c506e9324572788c4c1e'}]; written: 2026-11-12T23:41:31 vs 2026-11-12T01:25:34 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_b.svg | differs | svg (date): 7 of 129 lines differ, first at line 43 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_c.png | differs | png pixels |
| 5_figures/sfig_samap/panels/sfig_samap_panel_c.provenance.json | differs | json: git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; host: holy8a25505.rc.fas.harvard.edu vs holy7c10412.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...edeab02db39347443335007a'}] vs [{'path': '<RUN>/4...5e45c506e9324572788c4c1e'}]; written: 2026-11-12T23:41:57 vs 2026-11-12T01:25:53 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_c.svg | differs | svg (date): 12 of 184 lines differ, first at line 43 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_c_600dpi.png | differs | png pixels |
| 5_figures/sfig_samap/panels/sfig_samap_panel_d.provenance.json | differs | json: git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; host: holy8a25505.rc.fas.harvard.edu vs holy7c10412.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...31e71e6d3ed791be721bea62'}] vs [{'path': '<RUN>/4...31e71e6d3ed791be721bea62'}]; written: 2026-11-12T23:42:21 vs 2026-11-12T01:26:10 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_d.svg | differs | svg (date): 13 of 256 lines differ, first at line 48 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_e.png | differs | png pixels |
| 5_figures/sfig_samap/panels/sfig_samap_panel_e.provenance.json | differs | json: git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; host: holy8a25505.rc.fas.harvard.edu vs holy7c10412.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...edeab02db39347443335007a'}] vs [{'path': '<RUN>/4...5e45c506e9324572788c4c1e'}]; spearman_rho/all human cells: 0.7623 vs 0.6659; spearman_rho/decidual cells: 0.7341 vs 0.7273; spearman_rho/fibroblasts: 0.5804 vs 0.3074; written: 2026-11-12T23:42:52 vs 2026-11-12T01:26:40 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_e.svg | differs | svg (date): 51 of 379 lines differ, first at line 43 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_e_600dpi.png | differs | png pixels |
| 5_figures/sfig_samap/panels/sfig_samap_panel_f.provenance.json | differs | json: git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; host: holy8a25505.rc.fas.harvard.edu vs holy7c10412.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...69a0057fcc45e2585054e83a'}] vs [{'path': '<RUN>/4...69a0057fcc45e2585054e83a'}]; written: 2026-11-12T23:43:22 vs 2026-11-12T01:27:01 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_f.svg | differs | svg (date): 11 of 142 lines differ, first at line 34 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_g.png | differs | png pixels |
| 5_figures/sfig_samap/panels/sfig_samap_panel_g.provenance.json | differs | json: comparison: [{'design': 'unifo...ved': 0, 'shuffles': 1000}] vs [{'design': 'unifo...ved': 0, 'shuffles': 1000}]; git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; host: holy8a25505.rc.fas.harvard.edu vs holy7c10412.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...edeab02db39347443335007a'}] vs [{'path': '<RUN>/4...5e45c506e9324572788c4c1e'}]; permutation_observed/rho_all: 0.7623 vs 0.6659; permutation_observed/rho_decidual: 0.7341 vs 0.7273; permutation_observed/rho_fibroblast: 0.5804 vs 0.3074; written: 2026-11-12T23:43:44 vs 2026-11-12T01:27:17 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_g.svg | differs | svg (date): 74 of 549 lines differ, first at line 43 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_g_600dpi.png | differs | png pixels |
| 5_figures/sfig_samap/sfig_samap.png | differs | png pixels |
| 5_figures/sfig_samap/sfig_samap.provenance.json | differs | json: git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; host: holy8a25505.rc.fas.harvard.edu vs holy7c10412.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/5...d513a38407a868bbbbb4c5ba'}] vs [{'path': '<RUN>/5...b288ac40d58cb8065e1df9c6'}]; panels/a/git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; panels/a/written: 2026-11-12T23:41:09 vs 2026-11-12T01:25:12; panels/b/git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; panels/b/written: 2026-11-12T23:41:31 vs 2026-11-12T01:25:34; panels/c/git/commit: 33a8f5ca1731f22a20396449beefcac66b536914 vs 19231933290a5f0856c20636307bd6ab7a1fe57c; ... 17 more |
| 5_figures/sfig_samap/sfig_samap.svg | differs | svg (date): 173 of 1,754 lines differ, first at line 43 |
| 5_figures/sfig_samap/sfig_samap_600dpi.png | differs | png pixels |
| compare_with_record_20261112_012939.csv | only in other |  |
| notebooks/1_inputs.ipynb | differs | notebook code: 1 of 5 code cells differ |
| notebooks/2_samap.ipynb | same content | notebook: same 8 code cells (outputs and run times not compared) |
| notebooks/3_transfer.ipynb | differs | notebook code: 1 of 5 code cells differ |
| notebooks/4_statistics.ipynb | differs | notebook code: 1 of 13 code cells differ |
| notebooks/5_figures/fig6_umap_transferred_cell_type.ipynb | same content | notebook: same 9 code cells (outputs and run times not compared) |
| notebooks/5_figures/fig6c_umap_transferred_pseudotime.ipynb | same content | notebook: same 9 code cells (outputs and run times not compared) |
| notebooks/5_figures/fig6e_sections_transferred_pseudotime.ipynb | same content | notebook: same 9 code cells (outputs and run times not compared) |
| notebooks/5_figures/sfig_samap_a_umap_mouse_pseudotime.ipynb | same content | notebook: same 9 code cells (outputs and run times not compared) |
| notebooks/5_figures/sfig_samap_b_umap_transferred_pseudotime.ipynb | same content | notebook: same 9 code cells (outputs and run times not compared) |
| notebooks/5_figures/sfig_samap_c_umap_human_pseudotime.ipynb | same content | notebook: same 9 code cells (outputs and run times not compared) |
| notebooks/5_figures/sfig_samap_compose.ipynb | same content | notebook: same 8 code cells (outputs and run times not compared) |
| notebooks/5_figures/sfig_samap_d_edge_share.ipynb | same content | notebook: same 9 code cells (outputs and run times not compared) |
| notebooks/5_figures/sfig_samap_e_transferred_vs_human.ipynb | same content | notebook: same 10 code cells (outputs and run times not compared) |
| notebooks/5_figures/sfig_samap_f_sections_shuffled.ipynb | same content | notebook: same 10 code cells (outputs and run times not compared) |
| notebooks/5_figures/sfig_samap_g_permutation_statistics.ipynb | same content | notebook: same 10 code cells (outputs and run times not compared) |

## Manifest

samap/record/manifests/2026_11_12_xmens_samap_final.sha256 (sha256 ef4410acd828d5607e59643dff4ff75aa395a6b6320da6d48d17ea16a1bdc2ff)
