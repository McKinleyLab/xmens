# record/: the SAMap run of record

Unedited copies of the scripts, job files, logs and small run outputs that made the SAMap results of the paper:
the fit behind Fig. 6C, 6E, the cell-identity panel and Supplementary Fig. S30, and the deposited
`final_objects/Samap_Final.h5ad`. Nothing here is run by the pipeline; the step scripts one folder up reimplement it,
and `compare_with_record` checks their outputs against the files listed in `OUTPUTS.csv` and the numbers in
`numbers.csv`.

| folder | what made the run of record | when |
|---|---|---|
| `0_homology/` | the homology maps SAMap reads: `phmmer_hs_vs_mm.sh` and `phmmer_mm_vs_hs.sh` (32-task arrays; HMMER 3.4 phmmer, `-E 1e-5`, 4 threads) searched each proteome's 32 chunks against the other proteome (Ensembl GRCh38 and GRCm39 pep.all); `phmmer_to_samap_symbols.py` turned the hits into `homology/hsmm_hmmer/hsmm/{hs_to_mm,mm_to_hs}.txt` (the `.bak` is the version that ran; the copy without it only gained path arguments on 10 Sep). The conversion job of record, 8045975, failed on a wrong script path and the run that wrote the maps left no log; the converter re-run on the hit tables of record gives byte-identical maps (10 Sep, `homology_rebuild_check.txt`; again 6 Nov). The chunks were a round-robin split of each proteome (chunk i holds proteins i, i + 32, i + 64, ..., counting from 1; checked 6 Nov). `phmmer_tables_sha256.csv` checksums the 64 hit tables; `19_homology.smk` is xmens_final's Snakemake description of the step, for reference only | 24 Apr 2026 |
| `1_inputs/` | `job_prep_raw.sbatch` running `prep_full_input.py --scale raw` on the netscratch `2026_08_28` fib/dec objects (not the `2026_09_15` release that `21_samap_raw.smk` names) | 1 Sep 2026, job 43728073 |
| `2_samap/` | `job_arms_raw.sbatch` running `run_arm.py`, arm `vanilla_raw`: raw counts, `min_expression` 1, `NUMITERS` 3, homology maps `hsmm_hmmer`, 8 cores, 128 GB. SAM.run 1.4 min (human) and 0.8 min (mouse), SAMap init 10.9 min, SAMap run 33.1 min; mean alignment 0.660 | 1 Sep 2026, job 43728765 (job 43728340 failed on a scipy.sparse error first) |
| `3_transfer/` | label transfer `label_stain_full.py` (job 43744446, 1 Sep); joint export `make_transferred_h5ad.py` (job 43987051, 2 Sep); pseudotime transfer of record made by hand on 24 Oct with no log, rebuilt identically on 3 Nov by `samap_transfer_pt.py` and `samap_attach_pt.py` (`transfer_rebuild_2026_11_03_README.md`); obs of the deposited object cut to six columns by `slim_samap_obs.py` (2 Nov) | 1 Sep - 2 Nov 2026 |
| `4_5_figures/` | `29_fig_samap.smk`, its scripts and notebooks, and `samap_figures_final_objects.yaml`, at xmens_final `59e00dd`; the final figure set was drawn with them on 3 Nov in `tmp/2026_11_03_fig6c_pipeline/repo` and checked against the 2 Nov render (`final_check_2026_11_03_README.md`). `results_final/` is a copy of that set (`samap_figures/`: tables, graph block, edge shares, permutation control; `samap_main/`: Fig. 6C, 6E and the cell-identity panel; `samap_supp/`: Supplementary Fig. S30), which `compare_with_record` compares runs with | 3 Nov 2026 |
| `config/` | SAMap environment spec and lock (`workflow/envs/samap`), and the package lists of `samap_harmony_uv` and `rsc_uv` as installed on 6 Nov 2026 | |
| `methods_check/` | numbers the Methods report that the SAMap scripts of record did not compute: the mean total cross-species edge weight per compartment (12.62-13.53 of at most 20), SAMap's alignment scores between compartments (0.58 and 0.58 matched, 0.08 and 0.12 off-diagonal) and the nuclei without a cross-species edge (24 human, 1 mouse). The Methods took them from the methods-numbers check of 29 Oct (`pipeline_numbers_samap_rows.csv`, six digits); `graph_numbers_full.py` (with its job file, log and output) recomputes them independently at full precision from the joint graph of the fit of record with the final labels, and gives the same values. The score is SAMap's formula (`samap.analysis.scores._compute_csim`), checked against SAMap's own `mapping_celltype.csv` on two graphs (the record fit with its fit-time labels, and the fresh fit of 6 Nov): equal within 1.1e-16. `ensembl_release_check.txt`: the proteomes are byte-identical to the decompressed Ensembl release-113 files | 7 Nov 2026 |

- `SOURCES.csv`: every file here, its source (a path, or `xmens_final@<commit>:<path>`) and sha256; each copy was
  checked equal to its source.
- `OUTPUTS.csv`: the large inputs and outputs of the run of record, checksummed in place (homology maps, the
  `2026_08_28` source objects, `inputs_raw/`, `samap.pkl`, the joint export, the transfer table and object, and the
  deposited objects).
- `numbers.csv`: every number the figures and methods cite (fit settings and metrics, transfer counts, label
  agreement, 6C and 6E counts, edge shares, the permutation statistics), each with the file it was read from.
- `LOGS_AND_OUTPUTS.md`: where each step's logs and outputs are.

The fit sets no seed of its own. The bundled SAM seeds its PCA, neighbor search and run at 0 by default, yet a fresh
fit on identical inputs differs from the record (6 Nov 2026: 72% of the record's graph edges; see the pipeline's
README). Where the differences arise has not been pinned down; the approximate neighbor index (hnswlib) is built
with several threads.
