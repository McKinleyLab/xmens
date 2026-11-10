# record/: the MI run of record (copy of record, 2026-11-06)

Unedited copies of the scripts, settings, logs and small outputs that made the pseudotime MI of the paper: the MI
tables (`xmens_final/results/mi_dpt/pt_harmonized_{mouse,human}_nall/mi_all.csv`), the Methods overlap and Fig. S29C-D
sweep (`overlap_of_record/`), the Fig. 5F gene selection, and the supplementary table approved on 3 Nov 2026. Nothing
here is run; the step scripts one folder up reimplement it, and `compare_with_record` checks their outputs against the
files listed in `OUTPUTS.tsv`.

| Step | Files here | Ran | Output |
|---|---|---|---|
| 1-3 MI, mouse | `1-3_mi_pseudotime/`: `mi_pseudotime.py` (prep, worker, merge), `mouse_worker.sbatch` (the array of record), `mouse_manifest.json`, `mouse_worker_task0.out` | workers 22 Oct 18:07 - 23 Oct 08:37, 2026 (SLURM array 54745394, 200 tasks of 3.7-6.1 h); prep files and `mi_all.csv` dated 26 Oct, when Snakemake adopted the arm | `results/mi_dpt/pt_harmonized_mouse_nall/` |
| 1-3 MI, human | `build_mi_record_input.py` (puts the 30 Oct refit pseudotime into the 22 Oct object), `mi_pseudotime.py`, `22_mi_pseudotime.smk` (rules `mi_record_input`, `mi_prep`, `mi_worker` x 1,500, `mi_merge`), `human_build_meta.json`, `human_manifest.json`, `human_prep.log`, `human_worker_000.log`, `human_merge.log` | 31 Oct 21:45 - 1 Nov 01:49, 2026, two Snakemake runs (input job 56782476, prep 56782734); 439 workers failed within a minute on two nodes and ran again in the second | `results/mi_dpt/pt_harmonized_human_nall/` |
| 4 Overlap | `4_overlap/`: `mi_overlap.py` (rule `mi_overlap`), `overlap.log`, `methods_overlap.json` | 1 Nov 2026, 01:51 | `results/mi_dpt/overlap_of_record/` |
| 5 Table | `5_table/`: `build_human_mouse_mi_table.py` and the approved workbook `SuppTable_pseudotime_MI.xlsx` | 3 Nov 2026, on a Mac, from the two `mi_all.csv` and `ortholog_pairs_one2one.csv` (sha256 in its README sheet) | Drive `supplementary_tables/human_mouse_mi_final_pseudotime.xlsx` (identical; session artifact version 5d2621f8) |
| 6 Fig. 5F genes | `6_fig5f/`: `crossspecies_gene_selection.py` (rule `mi_heatmap_selection`), `heatmap_selection.log`, `selection_summary.json`, `conserved_profiles_selected_genes.csv`; `fig5f_panel_genes.csv`, the 28 genes drawn | 1 Nov 2026, 01:55 | `results/mi_dpt/overlap_of_record/crossspecies_heatmap_selection/` |

Paths in the table are relative to `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final`.

- Settings: `config/config.yaml`, blocks `seed` (line 1), `mi` (76-96), `pseudotime_of_record` (166-181) and
  `mi_pseudotime` (323-377). The copy is the working file of 3 Nov; its only change from the commit (d0983ff) is an added
  `ivd:` block.
- `SOURCES.tsv`: every file here, its source and sha256; each copy was checked with `cmp`. The two `5_table/` files come
  from this analysis session's artifacts, not from a cluster path.
- `OUTPUTS.tsv`: the large inputs and outputs of the run of record, checksummed in place on 6 Nov 2026 (the per-chunk
  results as one sha256 over the files concatenated in name order).
- `LOGS_AND_OUTPUTS.md`: where each step's logs and outputs are.

Code versions. The scripts here are xmens_final commit 39c83cf (31 Oct 2026). Their worker and merge code are identical
to commit 00724a4 (15 Sep). The mouse workers ran on 22 Oct, perhaps with the "readability" version of 20 Oct (e5d1c21),
whose worker and merge differ only in variable names: same seeds, same shared shuffles, float32 null, same p and merge.
The mouse prep ran an earlier prep (its manifest lacks `min_counts_source` and `n_cells_dropped_missing_pt`) on
`results/mi_dpt/pt_harmonized_mouse/source_subset.h5ad`, which no longer exists; step 1 here reproduces its matrices
exactly (README.md one folder up).

Known breaks, if anyone reruns these copies:
- The `*.smk` file is a Snakemake rule file; `mi_pseudotime.py` and the other scripts are run by it with the settings of
  `config.yaml`, and read and write inside `xmens_final/`.
- `results/mi_dpt/pt_harmonized_human_nall/worker.sbatch` is a 22 Oct file for an earlier 200-chunk human run; the human
  arm of record ran 1,500 chunks under Snakemake (rule `mi_worker`), so it is not copied here.
- `mi_overlap.py` also needs the two mouse control arms (`pt_harmonized_unind_mouse_nall`,
  `pt_harmonized_indmatch_mouse_nall`), which this folder does not compute.
- `build_human_mouse_mi_table.py` takes `--mi-mouse`, `--mi-human`, `--pairs`, `--record-overlap` and `--out`; it was
  run with the cluster paths given as labels (`--mi-mouse-label` and so on) on local copies of the same files.
