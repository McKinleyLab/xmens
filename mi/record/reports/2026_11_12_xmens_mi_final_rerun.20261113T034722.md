# mi: 2026_11_12_xmens_mi_final_rerun

Written 2026-11-13 03:47:22 by record/runs_report.py

OUT: /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_xmens_mi_final_rerun

## Steps (sacct)

- 1_inputs 60199204: COMPLETED (exit 0:0), 0:03:51, holy8a28512, 2026-11-12 23:12 to 2026-11-12 23:16
- 2_mi 60199207: 1,700 tasks: 1,700 COMPLETED (exit 0:0), 0:30:46-4:25:22 per task, 37 nodes, 2026-11-12 23:16 to 2026-11-13 03:43
- 3_combine 60199209: COMPLETED (exit 0:0), 0:00:33, holy8a27408, 2026-11-13 03:43 to 2026-11-13 03:43
- 4_overlap 60199230: COMPLETED (exit 0:0), 0:00:14, holy8a27408, 2026-11-13 03:44 to 2026-11-13 03:44
- 5_table 60199232: COMPLETED (exit 0:0), 0:00:22, holy8a27408, 2026-11-13 03:44 to 2026-11-13 03:44
- 6_fig5f_genes 60199234: COMPLETED (exit 0:0), 0:02:38, holy8a25503, 2026-11-13 03:44 to 2026-11-13 03:47
  - 2_mi: states {'COMPLETED': 1700}, exit codes {'0:0': 1700}, nodes holy7c10411, holy7c10412, holy7c10501, holy7c10502, holy7c10503, holy8a25503, holy8a25504, holy8a25505, holy8a25506, holy8a25603, holy8a25604, holy8a25605, holy8a25606, holy8a27111, holy8a27112, holy8a27207, holy8a27208, holy8a27209, holy8a27210, holy8a27211 ...

## Outputs

1,728 files, 292.7 MB; sha256 of each in mi/record/manifests/2026_11_12_xmens_mi_final_rerun.sha256 (sha256 f2c9643ddebc77f9...)

## Numbers

- N: 9532
- K: 1133
- n: 1215
- k: 375
- log10_p: -82.60218599263335

## Compared with /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_xmens_mi_final

1,728 files: 1,715 identical bytes; 5 same content (3 .json, 2 .gz; time stamps or paths only); 7 logs; 1 differ; 0 in one folder only

Files not identical in bytes, by folder:

- fig5f: 3 same content
- logs: 7 log
- overlap: 1 same content
- table: 1 differs, 1 same content

| file | result | detail |
|---|---|---|
| fig5f/percell_human_pool.tsv.gz | same content | gzip (header time stamp): identical once each run's folder path is replaced by <RUN> |
| fig5f/percell_mouse_pool.tsv.gz | same content | gzip (header time stamp): identical once each run's folder path is replaced by <RUN> |
| fig5f/selection_summary.json | same content | json: identical once each run's folder path is replaced by <RUN> |
| logs/1_inputs_human.log | log | logs record times and hosts; not compared |
| logs/1_inputs_mouse.log | log | log in one folder only |
| logs/3_combine_human.log | log | logs record times and hosts; not compared |
| logs/3_combine_mouse.log | log | logs record times and hosts; not compared |
| logs/4_overlap.log | log | logs record times and hosts; not compared |
| logs/5_table.log | log | logs record times and hosts; not compared |
| logs/6_fig5f_genes.log | log | logs record times and hosts; not compared |
| overlap/methods_overlap.json | same content | json: identical once each run's folder path is replaced by <RUN> |
| table/SuppTable_pseudotime_MI.xlsx | differs | xlsx sheets differ: README |
| table/SuppTable_pseudotime_MI_provenance.json | same content | json: differs only in built |

## Manifest

mi/record/manifests/2026_11_12_xmens_mi_final_rerun.sha256 (sha256 f2c9643ddebc77f99ad95785229be9bf569816791be84c46a52237a3a69a3649)
