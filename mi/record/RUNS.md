
## 2026-11-12 — MI, human side on the final human object (eFib1-root pseudotime)
Purpose:  human MI genes for the final human pseudotime; mouse side unchanged, reused from the November run
Code:     mi/ @ 33a8f5c
Inputs:   H5AD_HUMAN /n/lab_storage/eddy_lab/Lab/mckinley/data/20261112_Human_fibdec_final_nick_pdt.h5ad
                     sha256 c7d4bacb40f6472000c372cd005c6712a9c3f0e112ae5b23fb7df04b31dea82b
          H5AD_MOUSE /n/lab_storage/eddy_lab/Lab/mckinley/data/20261110_mouse_only_eFib_dec_final_nick_pdt.h5ad
                     sha256 e31372ba5fdd706ad620e8259e34216bc96b8f3b6394e4d2eedc0baeec8760dd
          mouse inputs and the 200 mouse chunks copied from
          /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_11_xmens_mi_november
OUT:      /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_xmens_mi_final
Jobs:     1_inputs 60180911 (human only); 2_mi 60192292 (array 200-1699, the human chunks); 3_combine 60192309;
          4_overlap 60192313; 5_table 60192330; 6_fig5f_genes 60192345 (mi/record/jobs.tsv)
Step 1:   human 44,949 nuclei x 15,354 genes, as in November; genes.txt and expression.npz byte-identical to
          November (cmp); pseudotime.npy differs, Spearman 0.7909 vs November (the root change)
<!-- report 2026_11_12_xmens_mi_final (record/runs_report.py fills in this block; written 2026-11-13 03:47) -->
Ran:      1_inputs 60180911: COMPLETED (exit 0:0), 0:03:43, holy7c10501, 2026-11-12 21:53 to 2026-11-12 21:57
          2_mi 60192292: 1,500 tasks: 1,500 COMPLETED (exit 0:0), 0:31:20-1:05:20 per task, 136 nodes, 2026-11-12 22:34
          to 2026-11-12 23:54
          3_combine 60192309: COMPLETED (exit 0:0), 0:00:29, holy7c10502, 2026-11-12 23:55 to 2026-11-12 23:56
          4_overlap 60192313: COMPLETED (exit 0:0), 0:00:12, holy7c10502, 2026-11-12 23:56 to 2026-11-12 23:56
          5_table 60192330: COMPLETED (exit 0:0), 0:00:31, holy8a27408, 2026-11-12 23:56 to 2026-11-12 23:57
          6_fig5f_genes 60192345: COMPLETED (exit 0:0), 0:01:59, holygpu8a26504, 2026-11-12 23:56 to 2026-11-12 23:58
Outputs:  1,727 files, 292.7 MB; sha256 of each in mi/record/manifests/2026_11_12_xmens_mi_final.sha256 (sha256
          49485cdcb75747df...)
Compared: with 2026_11_11_xmens_mi_november (reference; differences expected and not judged), file by file: 1,735 files:
          208 identical bytes; 8 logs; 1,513 differ; 6 in one folder only
Rerun:    2026_11_12_xmens_mi_final_rerun: outputs DIFFER (see that entry)
Numbers:  N 9532; K 1133; n 1215; k 375; log10_p -82.60218599263335
Code now: mi/ as at 33a8f5c
Report:   mi/record/reports/2026_11_12_xmens_mi_final.20261113T034722.md
Status:   complete; NOT reproduced by its rerun
<!-- end report 2026_11_12_xmens_mi_final -->

## 2026-11-12 — MI rerun, both species, of the final run (2026_11_12_xmens_mi_final)
Purpose:  every MI step again from scratch, mouse and human, to show the run of record reproduces
Code:     mi/ @ 33a8f5c
Inputs:   as the final run (H5AD_HUMAN and H5AD_MOUSE, sha256 pinned in config.sh)
OUT:      /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_xmens_mi_final_rerun
Jobs:     mi/record/jobs.tsv (submitted with --nice=10000, behind the run of record)
<!-- report 2026_11_12_xmens_mi_final_rerun (record/runs_report.py fills in this block; written 2026-11-13 03:47) -->
Ran:      1_inputs 60199204: COMPLETED (exit 0:0), 0:03:51, holy8a28512, 2026-11-12 23:12 to 2026-11-12 23:16
          2_mi 60199207: 1,700 tasks: 1,700 COMPLETED (exit 0:0), 0:30:46-4:25:22 per task, 37 nodes, 2026-11-12 23:16
          to 2026-11-13 03:43
          3_combine 60199209: COMPLETED (exit 0:0), 0:00:33, holy8a27408, 2026-11-13 03:43 to 2026-11-13 03:43
          4_overlap 60199230: COMPLETED (exit 0:0), 0:00:14, holy8a27408, 2026-11-13 03:44 to 2026-11-13 03:44
          5_table 60199232: COMPLETED (exit 0:0), 0:00:22, holy8a27408, 2026-11-13 03:44 to 2026-11-13 03:44
          6_fig5f_genes 60199234: COMPLETED (exit 0:0), 0:02:38, holy8a25503, 2026-11-13 03:44 to 2026-11-13 03:47
Outputs:  1,728 files, 292.7 MB; sha256 of each in mi/record/manifests/2026_11_12_xmens_mi_final_rerun.sha256 (sha256
          f2c9643ddebc77f9...)
Checks:   with 2026_11_12_xmens_mi_final, file by file: 1,728 files: 1,715 identical bytes; 5 same content (3 .json, 2
          .gz; time stamps or paths only); 7 logs; 1 differ; 0 in one folder only
Numbers:  N 9532; K 1133; n 1215; k 375; log10_p -82.60218599263335
Code now: mi/ as at 33a8f5c
Report:   mi/record/reports/2026_11_12_xmens_mi_final_rerun.20261113T034722.md
Status:   complete; DIFFERS from the run of record (see the report)
<!-- end report 2026_11_12_xmens_mi_final_rerun -->

## 2026-11-13 — MI steps 4-6, ortholog pairs matched by Ensembl gene ID
Purpose:  the MI overlap, supplementary table and Fig. 5F genes of record. The 12 Nov run named ortholog pairs by
          their Ensembl 116 symbols and joined them to the MI tables by name, which left out genes renamed since the
          objects' annotation (mouse Sars, now Sars1) and paired PALM2AKAP2 with Pakap instead of Pakap-1. Pairs are
          now matched to the objects' genes by Ensembl gene ID, by symbol where the ID changed. Supersedes steps 4-6
          of 2026_11_12_xmens_mi_final.
Code:     mi/ @ 4bbaf54 for 4_overlap and 6_fig5f_genes (submitted 10:04, two minutes after the commit); the 5_table
          rerun @ 0f407dd (short README; submitted 10:38:10, eight seconds after the commit). The run records no
          script sha256; each script on disk equals its committed version and was last modified before its job
          started.
Inputs:   mi_human.csv and mi_mouse.csv (step 3) and inputs/{human,mouse}/pseudotime.npy copied from
          2026_11_12_xmens_mi_final (identical, cmp)
          ORTHOLOGS (mi/config.sh), sha256 7c69c3a567f7e94dca9ef0a42be94da29ce196b9bb2bc8db0339339db654cdb7:
          /n/lab_storage/eddy_lab/users/nhilgert/human_ivd_repro/orthogroups/ensembl_human_mouse_components.tsv
          H5AD_HUMAN and H5AD_MOUSE (their gene IDs) as the 12 Nov run (sha256 c7d4bacb..., e31372ba...)
OUT:      /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_13_xmens_mi_final_ids
Jobs:     4_overlap 60288053; 5_table 60288059 (long README, superseded; its table moved to
          /n/lab_storage/eddy_lab/users/nhilgert/tmp/2026_11_13_mi_ids/table_long_readme); 6_fig5f_genes 60288061;
          5_table 60291933, the table of record (mi/record/jobs.tsv)
Step 4:   16,314 one-to-one components, 16,246 pairs. Ensembl ID not in the object: human 7 (3 matched by symbol),
          mouse 79 (12 by symbol). Object name differs from the Ensembl 116 symbol: human 147, mouse 627.
Table:    table/SuppTable_pseudotime_MI.xlsx, sha256 96fb0af44e6e57791aa82082c8b136fd706d5fae7d19f0af5d556250d0d71006;
          sent 13 Nov as SuppTable_pseudotime_MI_20261113.xlsx (the same file, renamed)
<!-- report 2026_11_13_xmens_mi_final_ids (record/runs_report.py fills in this block) -->
<!-- end report 2026_11_13_xmens_mi_final_ids -->
