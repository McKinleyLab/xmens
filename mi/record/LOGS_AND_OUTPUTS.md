# MI of record: where every script, log and output is (checked 2026-11-06)

X  = /n/lab_storage/eddy_lab/users/nhilgert/xmens_final
NS = /n/netscratch/eddy_lab/Everyone/nhilgert   (purged after 90 days)

## 1-3. Mouse MI arm, pt_harmonized_mouse_nall (22-23 Oct 2026)
- Scripts: X/workflow/scripts/mi/mi_pseudotime.py (prep, worker, merge); the array script
  X/results/mi_dpt/pt_harmonized_mouse_nall/worker.sbatch (`--array=0-199%100`, 32 GB, 14 h, `--perms 1001
  --subsample-cells 0`)
- Prep input: X/results/mi_dpt/pt_harmonized_mouse/source_subset.h5ad (manifest.json); the file no longer exists
- Workers: SLURM array 54745394, 22 Oct 18:07 - 23 Oct 08:37; 200 tasks completed, 3.7-6.1 h and 2.7 GB each (sacct);
  logs X/logs/mi_dpt/pt_harmonized_mouse_nall_<task>.out (two lines each)
- Merge: run when Snakemake adopted the arm on 26 Oct (mi_all.csv 14:35); X/logs/mi_dpt/pt_harmonized_mouse_nall_prep.log
  and _merge.log are empty
- Outputs: X/results/mi_dpt/pt_harmonized_mouse_nall/ expr.npz, genes.npy, pseudotime.npy, manifest.json,
  results/mi_000.csv-mi_199.csv, mi_all.csv

## 1-3. Human MI arm, pt_harmonized_human_nall (31 Oct - 1 Nov 2026)
- Rules: X/workflow/rules/22_mi_pseudotime.smk, mi_record_input -> mi_prep -> mi_worker (1,500) -> mi_merge; settings
  X/config/config.yaml (mi_pseudotime, pseudotime_of_record.human)
- Snakemake logs (X/.snakemake/log/):
  - 2026-10-31T214349.404936: submitted the input job (56782328), then stopped
  - 2026-10-31T214559.878451: input 56782476, prep 56782734, 1,061 of the 1,500 workers; 439 workers failed within a
    minute on nodes holy8a28511/12; ended 1 Nov 00:22 with "At least one job did not complete successfully"
  - 2026-11-01T002519.445262: the 439 remaining workers, the merge, mi_overlap and mi_heatmap_selection; ended 01:56
- Step logs: X/logs/mi_dpt/pt_harmonized_human_nall_prep.log, _worker_000.log to _worker_1499.log, _merge.log; SLURM logs
  under X/.snakemake/slurm_logs/rule_<rule>/ (Snakemake deletes those older than 10 days)
- Worker runtimes: at most 1.0 h each (median 0.76 h), 8 GB requested (sacct)
- Inputs: X/data/deliverables/2026_10_22_flagged_harmonized/Cevrim_Human_Only_Fib_Dec_Cells.h5ad and the pseudotime of
  record X/data/interim/human_pt_2026_10_30_run4_a40/cells.parquet (column pt_rank, sha256 956a07cd...)
- Outputs: X/results/mi_dpt/pt_harmonized_human_nall/ input.h5ad, build_meta.json, expr.npz, genes.npy, pseudotime.npy,
  manifest.json, results/mi_000.csv-mi_1499.csv, mi_all.csv
- Stale: X/results/mi_dpt/pt_harmonized_human_nall/worker.sbatch (22 Oct, 200 chunks) was not used for this arm

## 4. Overlap (1 Nov 2026, 01:51; Snakemake run 2026-11-01T002519)
- Script: X/workflow/scripts/mi/mi_overlap.py (rule mi_overlap); log X/logs/mi_dpt/overlap.log
- Inputs: the two arms above, the two mouse control arms (pt_harmonized_unind_mouse_nall, pt_harmonized_indmatch_mouse_nall,
  run 26 Oct), and /n/lab_storage/eddy_lab/users/nhilgert/human_ivd_repro/orthogroups/ensembl_human_mouse_components.tsv
- Outputs: X/results/mi_dpt/overlap_of_record/ ortholog_pairs_one2one.csv, methods_overlap.json,
  overlap_percentile_sweep.csv, and the tables not ported here (overlap_summary.csv, overlap_permutation_null.csv,
  overlap_cutoff_sweep.csv, shared_top_decile_genes.json, specificity_shared_universe.csv, overlap_report.json)

## 5. Supplementary table (3 Nov 2026)
- Script: build_human_mouse_mi_table.py, run on a Mac in this analysis session on copies of the two mi_all.csv and
  ortholog_pairs_one2one.csv (sha256 in the workbook's README sheet: 7351f3fe..., 6bcfc0a3..., 8e00db6a...)
- Output: the approved workbook, record/5_table/SuppTable_pseudotime_MI.xlsx (session artifact
  human_mouse_mi_final_pseudotime.xlsx, version 5d2621f8: MI in bits, results in bold, p written 1.95 x 10^-84)
- Drive (L/xmens_figures/supplementary_tables/, L = /n/lab_storage/eddy_lab/users/nhilgert):
  - human_mouse_mi_final_pseudotime.xlsx (3 Nov 11:10): the approved workbook, identical sha256 (ad1536da...)
  - SuppTable_pseudotime_MI.xlsx (2 Nov 23:51): an earlier build with the same data sheets, cell for cell; its overlap
    sheet also gives the thresholds in nats, and its checks sheet writes p as 2.0 x 10^-84. Its CSV exports are in
    csv/SuppTable_pseudotime_MI__*.csv and the script and provenance of that build in scripts_and_provenance/

## 6. Fig. 5F gene selection (1 Nov 2026, 01:55; Snakemake run 2026-11-01T002519)
- Script: X/workflow/scripts/mi/crossspecies_gene_selection.py (rule mi_heatmap_selection); log
  X/logs/mi_dpt/heatmap_selection.log
- Inputs: the two mi_all.csv, overlap_of_record/ortholog_pairs_one2one.csv, the mouse object
  NS/2026_10_22_flagged_harmonized/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad (obs pt_harmonized_rank) and the human arm's
  input.h5ad; X/resources/methods_numbers/fig5f_panel_genes.csv (the 28 genes drawn) for comparison
- Outputs: X/results/mi_dpt/overlap_of_record/crossspecies_heatmap_selection/ conserved_profiles_selected_genes.csv,
  conserved_profiles_pool_concordance.csv, percell_mouse_pool.tsv.gz, percell_human_pool.tsv.gz, selection_summary.json
