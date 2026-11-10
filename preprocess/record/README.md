# xmens_preprocess: FASTQ -> Cell Ranger -> CellBender -> QC (copy of record, 2026-11-06)

Unedited copies of the scripts that produced the February-regen QC'd objects
(/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_qc_regen_nick_andy_coords/qc_passed_{mouse,human}.h5ad).
SOURCES.tsv gives each file's original path, modification time and sha256 (each copy was checked with cmp).
Not runnable from here yet: paths inside the scripts still point to their original locations.

| Step | Files | Ran | Output |
|---|---|---|---|
| 0 Mouse reference | 0_reference/mus.sh (cellranger 7.2.0 mkref) | 2024-03-01 | reference/transgenic_mus |
| 1 Cell Ranger 7.2.0 count | 1_cellranger/: cellranger.sh, submit_cellranger.sh, submit_cellranger_human.sh; exact per-lane commands in feb_cmdlines/ | mouse 2026-02-18, human 2026-02-28 | decidualization/revisions/velocity/<lane>/cellranger/outs/ |
| 2 CellBender 0.2.2 | 2_cellbender/01_cellbender.smk (remove-background --total-droplets-included 50000 --cuda; env cellbender_022) | 2026-07-01 | decidualization/xmens_inputs_regen/<sample>/ |
| 3 Ingest + QC | 3_qc/: ingest.py, qc.py, qc_report.py, concat.py (xmens_final commit 59e00dd) and their rules (*.smk), run outside Snakemake by run_ingest_qc.py / run_ingest_qc.sbatch | 2026-11-06 | /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_qc_regen/ |
| 4 Positions | 4_positions/add_coords.py; Nick's coords from the Slide-tags fork (/n/lab_storage/eddy_lab/users/nhilgert/Slide-tags, separate git repo); Andy's from andy_raw_data/<lane>/Positions/ | 2026-11-06 | /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_qc_regen_nick_andy_coords/ |

Settings: config/config.yaml (qc: min_counts 500, n_mads 3, max_pct_mito 8; doublet: DoubletDetection 4.2, n_iters 40,
voter_thresh 0.5, p_thresh 1e-7, n_top_var_genes 3000) and config/samples.tsv.
Counts at every step, per lane: step_counts.tsv. Where the logs and outputs are: LOGS_AND_OUTPUTS.md.

How step 3 was launched (from /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_05_qc_verify_scripts_v2):

    OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_qc_regen
    J=$(sbatch --parsable --job-name=qc_regen --array=0-14 --export=ALL,INPUTS=regen,OUT=$OUT run_ingest_qc.sbatch)
    sbatch --dependency=afterok:$J --job-name=qc_regen_concat --export=ALL,INPUTS=regen,OUT=$OUT run_ingest_qc.sbatch

Step 1 re-run from the dataverse FASTQs: 1_cellranger/run_cellranger_count.sbatch (2026-11-06, SLURM array 58258549).
Step 4: add_coords.py was run with /n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python; its source and
destination folders are set inside the script.

Known breaks, to fix when detaching from Snakemake:
- The February scripts and feb_cmdlines use /n/eddy_lab/..., which no longer exists, and
  --fastqs=decidualization/lanes/<lane>/gex, now an empty folder (FASTQs: xmens_final/dataverse/fastqs/).
- cellranger.sh was last edited 2026-02-28 (human runs); for the mouse lanes, feb_cmdlines/ is the record.
- 01_cellbender.smk and the *.smk files are Snakemake rules; ingest.py, qc.py, qc_report.py and concat.py read a
  `snakemake` object, which run_ingest_qc.py builds.
- run_ingest_qc.py reads its config and scripts from /n/lab_storage/eddy_lab/users/nhilgert/xmens_final, not from here.
- run_ingest_qc.sbatch's final checks call three comparison scripts that are not copied here
  (compare_qc_outputs.py, compare_qc_vs_2025.py, check_qc_vs_submission.py).
- No script for adding GsD/GqD to the mouse FASTA/GTF, and no mkref log.
