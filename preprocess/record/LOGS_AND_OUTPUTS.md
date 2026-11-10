# Regen objects: where every script and log is (checked 2026-11-06)

L  = /n/lab_storage/eddy_lab/users/nhilgert
NS = /n/netscratch/eddy_lab/Everyone/nhilgert   (purged after 90 days)

## 0. Mouse reference (built 2024-03-01)
- Script: L/reference/mus.sh (cellranger 7.2.0 mkref)
- Inputs: L/reference/transgenic_mus_premkref/ (fasta/genome.fa, genes/genes.gtf, GsD_Berdeaux.fa, GqD_transgene.fa)
- Output: L/reference/transgenic_mus/ (reference.json records mkref 7.2.0)
- Logs: none found (mus.sh wrote mus_<jobid>.out/.err to its submit folder)
- Gap: no script found for adding the transgene contigs to genome.fa / genes.gtf
- Note: transgenic_mus/genes/ changed on 2026-02-17; it now holds an uncompressed genes.gtf dated 2024-03-01
- Human: L/reference/human/refdata-gex-GRCh38-2024-A (10x prebuilt)

## 1. FASTQ -> Cell Ranger 7.2.0 (mouse Feb 18, human Feb 28, 2026)
- Scripts: L/decidualization/revisions/velocity/cellranger.sh (per lane), submit_cellranger.sh (mouse),
  submit_cellranger_human.sh (human; contains rm -rf)
- cellranger.sh was last edited Feb 28 (human runs); each lane's exact command is in
  velocity/<lane>/cellranger/_cmdline (+ _invocation, _versions)
- Logs: velocity/<lane>/logs/<lane>_<jobid>.out/.err (2-6 jobs per lane);
  pipestance logs velocity/<lane>/cellranger/_log, _perf, cellranger.mri.tgz
- Recorded paths use /n/eddy_lab/..., which no longer exists; --fastqs points to
  decidualization/lanes/<lane>/gex, now an empty folder (FASTQs are in L/xmens_final/dataverse/fastqs/)
- Re-run (today): NS/2026_11_06_cellranger_count/run_cellranger_count.sbatch (twin L/tmp/2026_11_06_cellranger_count/)
  SLURM array 58258549; logs NS/2026_11_06_cellranger_count/logs/cellranger_count_58258549_<task>.out/.err

## 2. Cell Ranger -> CellBender 0.2.2 (Jul 1, 2026)
- Script: L/xmens_final/workflow/rules/01_cellbender.smk (defined only when ingest.rerun: true);
  config L/xmens_final/config/config.yaml lines 17-22; env L/.conda/envs/cellbender_022
- Logs: L/xmens_final/logs/00_cellbender_<sample>.log (15);
  L/decidualization/xmens_inputs_regen/<sample>/<sample>_cellbender.log (15, with command line) and _cellbender.pdf;
  Snakemake run log L/xmens_final/.snakemake/log/2026-07-01T122542.844634.snakemake.log
  (earlier attempts: 2026-06-30T230110.466922 and 2026-06-30T230608.783493)
- Gone: SLURM job logs .snakemake/slurm_logs/rule_cellbender/ (folder no longer exists)
- Stray: L/decidualization/xmens_inputs_regen/mck_1/ (Jun 30 test run, not used). Same cell barcodes as
  1_dpi_early_mck1/ (identical sha256); the filtered .h5 files differ

## 3. CellBender -> ingest + QC (Nov 6, 2026)
- Scripts (repo, commit 59e00dd): L/xmens_final/workflow/scripts/ingest_qc/ingest.py, qc.py, concat.py, qc_report.py;
  rules 03_ingest.smk, 04_qc.smk, 05_qc_report.smk, 06_concat.smk; config qc/doublet blocks in config.yaml
- Run outside Snakemake by: NS/2026_11_05_qc_verify_scripts_v2/run_ingest_qc.py + run_ingest_qc.sbatch
  (+ verify_qc.sh, compare_qc_outputs.py, compare_qc_vs_2025.py, check_qc_vs_submission.py);
  twin L/tmp/2026_11_05_qc_verify_scripts_v2/ (scripts only)
- Envs: L/.uv/envs/rsc_uv (ingest, concat), L/.conda/envs/cuda11.8_env (qc, DoubletDetection)
- Logs: NS/2026_11_06_qc_regen/logs/ (32), manifests/ (16 JSON: commit, envs, input sha256), specs/ (32),
  reports compare_20261106_063447.txt, vs_2025_*;
  SLURM outs NS/2026_11_05_qc_verify_scripts_v2/slurm-qc-rerun-58262416.out, -58262417.out and 15 task files
  (netscratch only; not in the twin)

## 4. Positions
- Nick's Slide-tags run (Mar 2026), git repo L/Slide-tags/: run_sample.sh, run_all_spatial.sh,
  batch_submit_sbcounts.sh, setup_fastq_dirs.sh, workflow/slidetag_pipe.sh, claude_edits/*.smk, samples/*.csv;
  output L/Slide-tags/coords/<lane>_coords.csv
  Logs: L/Slide-tags/data/<lane>/log/ (spatial_submit.log, 0run-4run/), data/<lane>/spatial/<index>/SBcounts/julia_run.log,
  L/Slide-tags/claude_status.log
- Andy's positions: L/../../Lab/mckinley_backup_2025-12-18/data/andy_raw_data/<lane>/Positions/ (no scripts or logs)
- Added to the objects by: L/tmp/2026_11_06_qc_regen_nick_andy_coords/add_coords.py (README_v2.md, tables)

## Final objects
- NS/2026_11_06_qc_regen_nick_andy_coords/qc_passed_{mouse,human}.h5ad (+ step_counts.tsv)

## Only on netscratch (save these)
1. NS/2026_11_06_qc_regen/: logs/, manifests/, specs/, reports (about 1 MB); raw/ and qc/ intermediates are 11 GB
2. NS/2026_11_05_qc_verify_scripts_v2/slurm-qc-rerun-*.out (17 files)
3. NS/2026_11_06_cellranger_count/ logs and pipestance files, once the run finishes
4. NS/2026_11_06_qc_regen_nick_andy_coords/ (the objects)
