# Logs and outputs of the SAMap run of record

Paths are on the harvard-eddy cluster. `FF` = `/n/lab_storage/eddy_lab/users/nhilgert/code/samap-scvi/full_fibdec`;
`NS` = `/n/netscratch/eddy_lab/Everyone/nhilgert` (purgeable). Checksums are in `OUTPUTS.csv`.

## Step 0: homology maps (24 Apr 2026)

`P` = `/n/lab_storage/eddy_lab/users/nhilgert/data/proteomes`; `SI` =
`/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/samap_integrate`. The scripts refer to both under
the retired `/n/eddy_lab` mount.

- Read: `P/fasta/Homo_sapiens.GRCh38.pep.all.fa` (123,845 proteins) and `P/fasta/Mus_musculus.GRCm39.pep.all.fa`
  (66,558), Ensembl, files dated Aug 2024 (release not recorded).
- Chunks: `P/chunks/{hs,mm}/<proteome>.{1..32}`, a round-robin split (chunk i holds proteins i, i + 32, i + 64, ..., counting from 1).
- phmmer (HMMER 3.4): `P/hm/hu_vs_mu/hu_<i>.tbl` and `P/hm/mu_vs_hu/mu_<i>.tbl` (12,833,563 and 12,780,006 hit lines),
  each with phmmer's full output beside it (`<prefix>_<i>.log`, about 1 GB, alignments included). `P/hm/*_merged.tbl`
  (the tables concatenated) and `P/hm/test/` are not read by anything. Task logs: `P/slurm_logs/` (copied to
  `0_homology/slurm_logs/`); about 50-55 min per task.
- Maps: `SI/homology/hsmm_hmmer/hsmm/hs_to_mm.txt` and `mm_to_hs.txt` (written 16:39-16:41, moved into `hsmm/` at
  17:15). The conversion job `SI/logs/phmmer_to_symbols_8045975.*` failed; the run that wrote the maps left no log.
  `SI/homology/hsmm_blast/` holds BLAST-based maps of the same day, not used by the fit of record.

## Step 1: SAMap inputs (1 Sep 2026)

- Read: `NS/2026_08_28_Cevrim_{Human,XMens}_Only_Fib_Dec_Cells.h5ad` (counts layer; 47,533 human and 43,197 mouse
  fib/dec nuclei).
- Wrote: `FF/inputs_raw/{hs,mm}_fib_dec_all_raw.h5ad` (raw integer counts in `.X`; human 36,601 genes).
- Log: `FF/logs/prep_raw_43728073.out` (copied to `1_inputs/`).

## Step 2: SAMap fit (1 Sep 2026)

- Wrote: `FF/runs_raw/vanilla_raw/samap.pkl` (10.7 GB), `eval.json`, `meta.json`, `mapping_*.csv`, `figures/`.
- Log: `FF/logs/samap_raw_43728765.out` (copied to `2_samap/`); the same job also fit the `harmony_raw` and
  `vanilla_raw_nofilter` arms, compared in `FF/runs_raw/summary.txt`.

## Step 3: transfers and the deposited object

- Label transfer (1 Sep): `FF/runs_raw/vanilla_raw/staining/{hs_from_mm_coarse,hs_from_mm_fine,mm_from_hs_coarse}.csv`;
  log `FF/logs/stain_all_43744446.out`.
- Joint export (2 Sep): `FF/runs_raw/vanilla_raw/samap_transferred.h5ad` (graph, UMAP and the label transfers in obs);
  log `FF/logs/th5_vanilla_raw_43987051.out`.
- Pseudotime transfer of record (24 Oct, no log): `FF/runs_raw/vanilla_raw/pt_transfer_harmonized/` and
  `FF/runs_raw/vanilla_raw/samap_transferred_pt_harmonized.h5ad` (sha256 `62fb75d7...`, `joint_sha256` in
  xmens_final `config.yaml`). Mouse source: `pt_harmonized` of
  `xmens_final/data/deliverables/2026_10_22_flagged_harmonized/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad`, rank-normalized.
  The rebuild by its rules (3 Nov) is in `/n/lab_storage/eddy_lab/users/nhilgert/tmp/2026_11_03_fig6c_pipeline/samap_raw_root/`.
- Deposited object (2 Nov): `NS/final_objects/Samap_Final.h5ad` (sha256 `ce15b9d4...`): the transfer object with obs
  cut to `species`, `sample`, `timepoint`, `cell_type`, `transferred_mouse_pseudotime`, `transferred_mouse_fibdec`;
  report `3_transfer/samap_final_slim_report.json`.

## Steps 4-5: statistics and figures (3 Nov 2026)

- Drawn by `snakemake samap_figures_all --configfile config/samap_figures_final_objects.yaml` in
  `/n/lab_storage/eddy_lab/users/nhilgert/tmp/2026_11_03_fig6c_pipeline/repo`:
  - `results/samap_figures/final/`: tables, graph block, edge shares, permutation control
    (`ptnull/ptnull_statistics_table.csv`, seed 20261022, 1,000 shuffles per design)
  - `results/figures/samap_main/final/`: Fig. 6C, 6E and the cell-identity panel, with provenance records
  - `results/figures/samap_supp/final/`: Supplementary Fig. S30, panels a-g and the assembled figure
- Inputs: the deposited fib/dec objects (`NS/final_objects/`, pinned in `4_5_figures/MANIFEST.json`) and
  `Samap_Final.h5ad`; the human pseudotime of record, `xmens_final/data/interim/human_pt_2026_10_30_run4_a40/cells.parquet`.
- On Drive: `xmens_figures/final_figures/Fig6_crossspecies_radial_patterning/`.
