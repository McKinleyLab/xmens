# samap: 2026_11_12_xmens_samap_final_rerun

Written 2026-11-13 03:47:22 by record/runs_report.py

OUT: /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_xmens_samap_final_rerun

## Steps (sacct)

- 1_inputs 60199554: COMPLETED (exit 0:0), 0:02:29, holygpu8a26504, 2026-11-12 23:09 to 2026-11-12 23:12
- 2_samap 60199556: COMPLETED (exit 0:0), 0:49:46, holy8a28511, 2026-11-12 23:12 to 2026-11-13 00:02
- 3_transfer 60199558: COMPLETED (exit 0:0), 0:02:38, holy8a27312, 2026-11-13 00:24 to 2026-11-13 00:27
- 4_statistics 60199560: COMPLETED (exit 0:0), 0:01:44, holy8a27312, 2026-11-13 00:27 to 2026-11-13 00:29
- 5_figures 60199581: COMPLETED (exit 0:0), 0:05:07, holy7c10501, 2026-11-13 00:29 to 2026-11-13 00:34

## Outputs

107 files, 15.0 GB; sha256 of each in samap/record/manifests/2026_11_12_xmens_samap_final_rerun.sha256 (sha256 0da8adf8dcfba6f9...)

## Numbers

- mean_alignment: 0.6731769968221357
- seed: 0
- hnswlib_threads: 1

## Compared with /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_xmens_samap_final

107 files: 51 identical bytes; 24 same content (15 .ipynb, 7 .pdf, 2 .json; time stamps or paths only); 1 not compared; 31 differ; 0 in one folder only

Files not identical in bytes, by folder:

- 1_inputs: 1 same content
- 2_samap: 1 same content, 1 not compared, 1 differs
- 3_transfer: 2 differs
- 4_statistics: 1 differs
- 4_statistics/graph_block: 1 differs
- 5_figures/fig6: 10 differs, 7 same content
- 5_figures/sfig_samap: 16 differs
- notebooks: 4 same content
- notebooks/5_figures: 11 same content

| file | result | detail |
|---|---|---|
| 1_inputs/inputs.json | same content | json: differs only in written |
| 2_samap/meta.json | same content | json: differs only in finished, host, started |
| 2_samap/samap.pkl | not compared | pickle; not loaded (its contents are compared through the h5ad outputs) |
| 2_samap/samap_joint.h5ad | differs | hdf5: 0 objects in one file only, 1 differ: uns/provenance |
| 3_transfer/Samap_Final.h5ad | differs | hdf5: 0 objects in one file only, 1 differ: uns/provenance |
| 3_transfer/summary.json | differs | json: samap_final/joint_sha256: 3ef81dd45dff9fbf54...bbbb51aec3d7f730730851302f7 vs b48e787da8cab11609...2c283748f0fab4c9ff8252e86eb |
| 4_statistics/graph_block/meta.json | differs | json: joint_sha256: 4197ad52e7146bbdda...81d4870541f951ded629e46fffe vs 7901c4992351a7be8b...624800982ae861d137363860e80 |
| 4_statistics/statistics.json | differs | json: graph_block/joint_sha256: 4197ad52e7146bbdda...81d4870541f951ded629e46fffe vs 7901c4992351a7be8b...624800982ae861d137363860e80; written: 2026-11-13T00:28:29 vs 2026-11-12T23:39:09 |
| 5_figures/fig6/20261027-Human-menstrual-transferred_pseudotime-spatial.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/20261027-Human-menstrual-transferred_pseudotime-spatial.svg | differs | svg (date): 28,377 of 42,562 lines differ, first at line 34 |
| 5_figures/fig6/20261027-Human-secretory-transferred_pseudotime-spatial.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/20261027-Human-secretory-transferred_pseudotime-spatial.svg | differs | svg (date): 36,043 of 54,077 lines differ, first at line 34 |
| 5_figures/fig6/20261027-spatial_key_transferred_pseudotime.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/20261027-spatial_key_transferred_pseudotime.svg | differs | svg (date): 5 of 86 lines differ, first at line 41 |
| 5_figures/fig6/fig6_umap_transferred_cell_type.provenance.json | differs | json: host: holy7c10501.rc.fas.harvard.edu vs holy8a25505.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...6cf3cf39a8678c1c042fd10d'}] vs [{'path': '<RUN>/4...6cf3cf39a8678c1c042fd10d'}]; written: 2026-11-13T00:34:29 vs 2026-11-12T23:46:05 |
| 5_figures/fig6/fig6c_umap_transferred_pseudotime.provenance.json | differs | json: host: holy7c10501.rc.fas.harvard.edu vs holy8a25505.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...d78495249e80edfe2cf8006b'}] vs [{'path': '<RUN>/4...d78495249e80edfe2cf8006b'}]; written: 2026-11-13T00:33:15 vs 2026-11-12T23:44:49 |
| 5_figures/fig6/fig6e_sections_transferred_pseudotime.provenance.json | differs | json: host: holy7c10501.rc.fas.harvard.edu vs holy8a25505.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...d78495249e80edfe2cf8006b'}] vs [{'path': '<RUN>/4...d78495249e80edfe2cf8006b'}]; written: 2026-11-13T00:33:38 vs 2026-11-12T23:45:20 |
| 5_figures/fig6/umap_dark_transferred_cell_type.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/umap_dark_transferred_cell_type.svg | differs | svg (date): 11 of 176 lines differ, first at line 41 |
| 5_figures/fig6/umap_dark_transferred_cell_type_humanonly.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/umap_dark_transferred_cell_type_humanonly.svg | differs | svg (date): 8 of 153 lines differ, first at line 43 |
| 5_figures/fig6/umap_dark_transferred_pseudotime.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/umap_dark_transferred_pseudotime.svg | differs | svg (date): 12 of 181 lines differ, first at line 41 |
| 5_figures/fig6/umap_dark_transferred_pseudotime_humanonly.pdf | same content | pdf: same apart from its dates and ID |
| 5_figures/fig6/umap_dark_transferred_pseudotime_humanonly.svg | differs | svg (date): 9 of 158 lines differ, first at line 43 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_a.provenance.json | differs | json: host: holy7c10501.rc.fas.harvard.edu vs holy8a25505.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...edeab02db39347443335007a'}] vs [{'path': '<RUN>/4...edeab02db39347443335007a'}]; written: 2026-11-13T00:30:24 vs 2026-11-12T23:41:09 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_a.svg | differs | svg (date): 5 of 125 lines differ, first at line 43 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_b.provenance.json | differs | json: host: holy7c10501.rc.fas.harvard.edu vs holy8a25505.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...edeab02db39347443335007a'}] vs [{'path': '<RUN>/4...edeab02db39347443335007a'}]; written: 2026-11-13T00:30:46 vs 2026-11-12T23:41:31 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_b.svg | differs | svg (date): 7 of 129 lines differ, first at line 43 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_c.provenance.json | differs | json: host: holy7c10501.rc.fas.harvard.edu vs holy8a25505.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...edeab02db39347443335007a'}] vs [{'path': '<RUN>/4...edeab02db39347443335007a'}]; written: 2026-11-13T00:31:07 vs 2026-11-12T23:41:57 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_c.svg | differs | svg (date): 12 of 184 lines differ, first at line 43 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_d.provenance.json | differs | json: host: holy7c10501.rc.fas.harvard.edu vs holy8a25505.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...31e71e6d3ed791be721bea62'}] vs [{'path': '<RUN>/4...31e71e6d3ed791be721bea62'}]; written: 2026-11-13T00:31:25 vs 2026-11-12T23:42:21 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_d.svg | differs | svg (date): 13 of 256 lines differ, first at line 48 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_e.provenance.json | differs | json: host: holy7c10501.rc.fas.harvard.edu vs holy8a25505.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...edeab02db39347443335007a'}] vs [{'path': '<RUN>/4...edeab02db39347443335007a'}]; written: 2026-11-13T00:31:49 vs 2026-11-12T23:42:52 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_e.svg | differs | svg (date): 39 of 379 lines differ, first at line 43 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_f.provenance.json | differs | json: host: holy7c10501.rc.fas.harvard.edu vs holy8a25505.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...69a0057fcc45e2585054e83a'}] vs [{'path': '<RUN>/4...69a0057fcc45e2585054e83a'}]; written: 2026-11-13T00:32:11 vs 2026-11-12T23:43:22 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_f.svg | differs | svg (date): 11 of 142 lines differ, first at line 34 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_g.provenance.json | differs | json: host: holy7c10501.rc.fas.harvard.edu vs holy8a25505.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/4...edeab02db39347443335007a'}] vs [{'path': '<RUN>/4...edeab02db39347443335007a'}]; written: 2026-11-13T00:32:29 vs 2026-11-12T23:43:44 |
| 5_figures/sfig_samap/panels/sfig_samap_panel_g.svg | differs | svg (date): 53 of 549 lines differ, first at line 43 |
| 5_figures/sfig_samap/sfig_samap.provenance.json | differs | json: host: holy7c10501.rc.fas.harvard.edu vs holy8a25505.rc.fas.harvard.edu; inputs: [{'path': '<RUN>/5...d2f5b399f4ee6c58e09dcfba'}] vs [{'path': '<RUN>/5...d513a38407a868bbbbb4c5ba'}]; panels/a/written: 2026-11-13T00:30:24 vs 2026-11-12T23:41:09; panels/b/written: 2026-11-13T00:30:46 vs 2026-11-12T23:41:31; panels/c/written: 2026-11-13T00:31:07 vs 2026-11-12T23:41:57; panels/d/written: 2026-11-13T00:31:25 vs 2026-11-12T23:42:21; panels/e/written: 2026-11-13T00:31:49 vs 2026-11-12T23:42:52; panels/f/written: 2026-11-13T00:32:11 vs 2026-11-12T23:43:22; ... 3 more |
| 5_figures/sfig_samap/sfig_samap.svg | differs | svg (date): 140 of 1,754 lines differ, first at line 43 |
| notebooks/1_inputs.ipynb | same content | notebook: same 5 code cells (outputs and run times not compared) |
| notebooks/2_samap.ipynb | same content | notebook: same 8 code cells (outputs and run times not compared) |
| notebooks/3_transfer.ipynb | same content | notebook: same 5 code cells (outputs and run times not compared) |
| notebooks/4_statistics.ipynb | same content | notebook: same 13 code cells (outputs and run times not compared) |
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

samap/record/manifests/2026_11_12_xmens_samap_final_rerun.sha256 (sha256 0da8adf8dcfba6f9a63f14ef8d1881da194f07dcdc30abd83d632bd3dbce5662)
