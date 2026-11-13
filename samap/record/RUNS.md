
## 2026-11-12 — SAMap steps 1-5 on the final objects (seeded fit)
Purpose:  SAMap object, transfers, statistics and figures for the final human (eFib1 root) and mouse objects
Code:     samap/ @ 33a8f5c (seeded fit since ca6f07d)
Inputs:   HUMAN_H5AD /n/lab_storage/eddy_lab/Lab/mckinley/data/20261112_Human_fibdec_final_nick_pdt.h5ad
                     sha256 c7d4bacb40f6472000c372cd005c6712a9c3f0e112ae5b23fb7df04b31dea82b
          MOUSE_H5AD /n/lab_storage/eddy_lab/Lab/mckinley/data/20261110_mouse_only_eFib_dec_final_nick_pdt.h5ad
                     sha256 e31372ba5fdd706ad620e8259e34216bc96b8f3b6394e4d2eedc0baeec8760dd
          MAPS /n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/samap_integrate/homology/hsmm_hmmer
               (homology maps of record, as in runs A and B; steps 0a/0b not run)
OUT:      /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_xmens_samap_final
Jobs:     1_inputs 60181191; 2_samap 60192327; 3_transfer 60192344; 4_statistics 60192359; 5_figures 60192360
          (samap/record/jobs.tsv)
Step 1:   mm_fib_dec_all_raw.h5ad byte-identical to run A (cmp); hs_fib_dec_all_raw.h5ad equal to run A in all 29
          datasets but uns/samap_input (source path); inputs.json differs only in the time and the human source
          (sha256 c7d4bacb...). Step-1 inputs hold no pseudotime: the step-2 fit should equal run A's.
<!-- report 2026_11_12_xmens_samap_final (record/runs_report.py fills in this block; written 2026-11-13 03:47) -->
Ran:      1_inputs 60181191: COMPLETED (exit 0:0), 0:02:13, holy8a25503, 2026-11-12 21:54 to 2026-11-12 21:56
          2_samap 60192327: COMPLETED (exit 0:0), 0:59:10, holy8a27311, 2026-11-12 22:35 to 2026-11-12 23:34
          3_transfer 60192344: COMPLETED (exit 0:0), 0:01:53, holy8a25505, 2026-11-12 23:35 to 2026-11-12 23:37
          4_statistics 60192359: COMPLETED (exit 0:0), 0:01:33, holy8a25505, 2026-11-12 23:38 to 2026-11-12 23:39
          5_figures 60192360: COMPLETED (exit 0:0), 0:05:43, holy8a25505, 2026-11-12 23:40 to 2026-11-12 23:46
Outputs:  107 files, 15.0 GB; sha256 of each in samap/record/manifests/2026_11_12_xmens_samap_final.sha256 (sha256
          ef4410acd828d560...)
Compared: with 2026_11_12_xmens_samap_seeded_A (reference; differences expected and not judged), file by file: 108
          files: 37 identical bytes; 19 same content (12 .ipynb, 7 .pdf; time stamps or paths only); 1 not compared; 50
          differ; 1 in one folder only
Rerun:    2026_11_12_xmens_samap_final_rerun: outputs DIFFER (see that entry)
Numbers:  mean_alignment 0.6731769968221357; seed 0; hnswlib_threads 1
Code now: differs from 33a8f5c in samap/README.md (records and logs aside)
Report:   samap/record/reports/2026_11_12_xmens_samap_final.20261113T034722.md
Status:   complete; NOT reproduced by its rerun
<!-- end report 2026_11_12_xmens_samap_final -->

## 2026-11-12 — SAMap rerun, steps 1-5, of the final run (2026_11_12_xmens_samap_final)
Purpose:  every SAMap step again from scratch (seeded fit), to show the run of record reproduces
Code:     samap/ @ 33a8f5c
Inputs:   as the final run (HUMAN_H5AD and MOUSE_H5AD, sha256 pinned in config.sh; MAPS as above)
OUT:      /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_12_xmens_samap_final_rerun
Jobs:     samap/record/jobs.tsv
<!-- report 2026_11_12_xmens_samap_final_rerun (record/runs_report.py fills in this block; written 2026-11-13 03:47) -->
Ran:      1_inputs 60199554: COMPLETED (exit 0:0), 0:02:29, holygpu8a26504, 2026-11-12 23:09 to 2026-11-12 23:12
          2_samap 60199556: COMPLETED (exit 0:0), 0:49:46, holy8a28511, 2026-11-12 23:12 to 2026-11-13 00:02
          3_transfer 60199558: COMPLETED (exit 0:0), 0:02:38, holy8a27312, 2026-11-13 00:24 to 2026-11-13 00:27
          4_statistics 60199560: COMPLETED (exit 0:0), 0:01:44, holy8a27312, 2026-11-13 00:27 to 2026-11-13 00:29
          5_figures 60199581: COMPLETED (exit 0:0), 0:05:07, holy7c10501, 2026-11-13 00:29 to 2026-11-13 00:34
Outputs:  107 files, 15.0 GB; sha256 of each in samap/record/manifests/2026_11_12_xmens_samap_final_rerun.sha256 (sha256
          0da8adf8dcfba6f9...)
Checks:   with 2026_11_12_xmens_samap_final, file by file: 107 files: 51 identical bytes; 24 same content (15 .ipynb, 7
          .pdf, 2 .json; time stamps or paths only); 1 not compared; 31 differ; 0 in one folder only
Numbers:  mean_alignment 0.6731769968221357; seed 0; hnswlib_threads 1
Code now: differs from 33a8f5c in samap/README.md (records and logs aside)
Report:   samap/record/reports/2026_11_12_xmens_samap_final_rerun.20261113T034722.md
Status:   complete; DIFFERS from the run of record (see the report)
<!-- end report 2026_11_12_xmens_samap_final_rerun -->

## 2026-11-13 — SAMap steps 0c and 2-5, homology maps under the objects' gene names (Ensembl-ID matching)
Purpose:  the SAMap object of record. The maps of record name genes by their Ensembl 113 symbols, so genes renamed
          since the objects' annotation (mouse Sept2, now Septin2) had no homology. Step 0c renames them in the maps
          to the objects' names, matched by Ensembl gene ID. Supersedes 2026_11_12_xmens_samap_final (same inputs and
          code otherwise), whose object is in LAB as 20261112_samap_final_gene_symbols.h5ad
          (sha256 7901c4992351a7be...).
Code:     samap/ @ b027750 (step 0c added; steps 2-5 as at 33a8f5c). Jobs submitted 10:05:16-20, two seconds after the
          commit; relabel.json records the sha256 of the script that ran (0e3c57e2...), the committed one. 5_figures
          ran after 0f407dd and records that commit, which changed mi/ only.
Inputs:   step 1 of 2026_11_12_xmens_samap_final (SAMAP_INPUTS=/n/netscratch/eddy_lab/Everyone/nhilgert/
          2026_11_12_xmens_samap_final/1_inputs)
          MAPS $OUT/0_homology_relabeled (step 0c), made from the maps of record
               /n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/samap_integrate/homology/hsmm_hmmer
               and the Ensembl 113 proteomes (sha256 pinned in config.sh)
OUT:      /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_13_xmens_samap_relabeled
Jobs:     0c_relabel 60288129; 2_samap 60288130; 3_transfer 60288131; 4_statistics 60288134; 5_figures 60288135
          (samap/record/jobs.tsv)
Step 0c:  human: 92 of 19,601 map names renamed, none left; mouse: 2,178 of 21,642 renamed (1,132 olfactory receptors,
          Or... in Ensembl 113, Olfr... in the object), 13 left (12 names held by two object genes, e.g. Pakap and
          Pakap-1; 1 name taken in the map). Lines changed: hs_to_mm.txt 530,224, mm_to_hs.txt 537,048; everything
          after the two gene names byte-identical.
Compared: with 2026_11_12_xmens_samap_final, from eval.json, transfers.csv and mapping_*.csv: mean_alignment 0.6732 ->
          0.6710; transferred cell type the same for 45,894 of 46,818 human nuclei (98.03%), the 924 others mostly
          uncertain calls (median vote share 0.60, against 1.00); transferred pseudotime Pearson r 0.9913; every human
          and mouse group keeps its best match in the cell type, subtype and timepoint mapping scores (largest
          change 0.0145).
LAB copy: 3_transfer/Samap_Final.h5ad copied to
          /n/lab_storage/eddy_lab/Lab/mckinley/data/20261113_samap_final_ensembl_ids.h5ad
          (sha256 a05a72788ee82ec90178d29b5116c0f4c340ab77aadcba954e0eb6ea2d187da3, identical); first posted in Slack
          as 20261112_samap_final_ensembl_ids.h5ad
<!-- report 2026_11_13_xmens_samap_relabeled (record/runs_report.py fills in this block) -->
<!-- end report 2026_11_13_xmens_samap_relabeled -->
