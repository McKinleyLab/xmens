# preprocess: FASTQs -> Cell Ranger -> CellBender -> QC'd objects

Five steps, one script each, no Snakemake. Paths live in `config.sh`; the 15 lanes are the rows of `samples.tsv`
(row n = array task n). Submit from this folder; SLURM output goes to `logs/`.

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/preprocess
export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_preprocess

sbatch 1_cellranger.sbatch     # Cell Ranger 7.2.0 count        -> $OUT/cellranger/<lane>/outs/
sbatch 2_cellbender.sbatch     # CellBender 0.2.2 (A40 GPU)      -> $OUT/cellbender/<sample>/
sbatch 3_qc.sbatch             # ingest + QC per sample          -> $OUT/raw/, $OUT/qc/
sbatch 4_concat.sbatch         # one object per species          -> $OUT/qc_passed_{mouse,human}.h5ad
sbatch 5_add_coords.sbatch     # adds Andy's positions to Nick's -> $OUT/with_coords/
```

Each step needs the previous one finished. To queue all five at once:

```bash
j1=$(sbatch --parsable 1_cellranger.sbatch)
j2=$(sbatch --parsable --dependency=afterok:$j1 2_cellbender.sbatch)
j3=$(sbatch --parsable --dependency=afterok:$j2 3_qc.sbatch)
j4=$(sbatch --parsable --dependency=afterok:$j3 4_concat.sbatch)
sbatch --dependency=afterok:$j4 5_add_coords.sbatch
```

- `bash 1_cellranger.sbatch --check` (also steps 2 and 3) prints what each task would run and whether its inputs and
  outputs exist, without running anything.
- One lane only: `sbatch --array=3 3_qc.sbatch` (row 3 of samples.tsv), or `SLURM_ARRAY_TASK_ID=3 bash 3_qc.sbatch`
  in an interactive session.
- Outputs are never overwritten: a task whose output exists is skipped, and files are written under a temporary
  name and renamed when complete.
- Each Python script also runs on its own; `python scripts/qc.py --help` lists its arguments.

## Starting from earlier results

Set these before submitting to skip steps 1-2 ({lane} and {sample} are filled in per task):

```bash
# the February 2026 Cell Ranger run (then start at step 2)
export CELLRANGER_RAW=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/velocity/{lane}/cellranger/outs/raw_feature_bc_matrix.h5
# the July 2026 CellBender run behind the regen objects (then start at step 3)
export CELLBENDER_H5=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/xmens_inputs_regen/{sample}/{sample}_cellbender_filtered.h5
```

## Settings

- Cell Ranger 7.2.0 `count`, 8 cores, 96 GB; mouse reference `reference/transgenic_mus` (mm10-2020-A plus GsD/GqD),
  human `refdata-gex-GRCh38-2024-A`.
- CellBender 0.2.2 `remove-background --total-droplets-included 50000 --cuda`, defaults otherwise.
- QC (`scripts/qc.py` defaults): at least 500 UMIs, at most 8% mitochondrial, 3 MADs on log10 total counts and
  genes; DoubletDetection 4.2 with 40 iterations, 3,000 genes, p_thresh 1e-7, voter_thresh 0.5; seed 0.
- Positions: Nick's March 2026 Slide-tags run (`Slide-tags/coords/<lane>_coords.csv`, attached in step 3) and Andy's
  (`andy_raw_data/<lane>/Positions/`, added in step 5).

## Differences from the run of record

`record/` holds the unedited scripts that made the November 2026 regen objects, with `SOURCES.tsv` (original paths
and sha256), `LOGS_AND_OUTPUTS.md` (where that run's logs are) and `step_counts.tsv` (counts at every step).
Compared with them:

- Snakemake inputs became command-line arguments; the computations are unchanged.
- `qc.py` stops if DoubletDetection is missing or fails, instead of skipping doublet removal with a warning.
- Cell Ranger adds `--no-bam` (no BAM; the matrices are unchanged) and `--disable-ui`.
- Step 1 avoids the nodes the February 2026 script excluded.
- Step 5 checks the objects' positions against Nick's CSVs before copying.

## Comparing with the run of record

```bash
sbatch compare_with_record.sbatch    # any time; compares the steps that have finished
```

The summary ends `logs/compare_with_record_<job ID>.out`; one row per sample goes to
`$OUT/compare_with_record_<date>_<time>.tsv`. It compares Cell Ranger with the February 2026 run (byte for byte),
CellBender with the July 2026 run (cells called; nuclei with identical counts), the QC numbers with the record's logs
(`2026_11_06_qc_regen/logs/01_qc_<sample>.log`), and the final objects and positions with the regen objects.

What reproduced in the checks of November 6, 2026:

- Cell Ranger: byte-identical matrices and metrics in the 11 lanes of the November re-run compared so far.
- Ingest, concat and positions: identical to the record (mck_1).
- QC thresholds: identical. DoubletDetection is not: on the same input, mck_1 had 70 doublets in the record and 68 in
  a repeat, all 68 among the 70.
- CellBender: two runs on mck_1 (June 30 and July 1, 2026, both on an A40 on holygpu7c0530 with the same
  environment) agree in every field of the output: counts, latent variables and training curve. The files differ
  only in HDF5 metadata bytes, so compare contents, not checksums.
