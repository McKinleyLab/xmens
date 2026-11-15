# preprocess: FASTQs -> Cell Ranger -> CellBender -> QC'd objects

One script per step. Paths are in `config.sh`; the 15 lanes are the rows of `samples.tsv` (row n = array task n).

## Run

```bash
cd preprocess && mkdir -p logs
export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_preprocess    # a new folder for each run

j1=$(sbatch --parsable 1_cellranger.sbatch)                              # Cell Ranger 7.2.0 count
j2=$(sbatch --parsable --dependency=afterok:$j1 2_cellbender.sbatch)     # CellBender 0.2.2 (A40 GPU)
j3=$(sbatch --parsable --dependency=afterok:$j2 3_qc.sbatch)             # ingest and QC per sample
j4=$(sbatch --parsable --dependency=afterok:$j3 4_concat.sbatch)         # one object per species
sbatch --dependency=afterok:$j4 5_add_coords.sbatch                      # adds the second set of positions
sbatch --dependency=afterok:$j3 6_qc_figures.sbatch                      # QC figures, one per sample
```

`bash 1_cellranger.sbatch --check` (also steps 2 and 3) prints what each task would run, without running anything.
One lane only: `sbatch --array=3 3_qc.sbatch`. A task whose output exists is skipped.

## Settings

- Cell Ranger 7.2.0 `count`; mouse reference `reference/transgenic_mus` (mm10-2020-A plus GsD/GqD), human
  `refdata-gex-GRCh38-2024-A`.
- CellBender 0.2.2 `remove-background --total-droplets-included 50000 --cuda`, defaults otherwise.
- QC (`scripts/qc.py` defaults): at least 500 UMIs, at most 8% mitochondrial, 3 MADs on log10 total counts and genes;
  DoubletDetection 4.2 with 40 iterations, 3,000 genes, p_thresh 1e-7, voter_thresh 0.5; seed 0. Step 3 runs on AMD
  nodes (`--constraint="genoa|turin"`), on which the doublet calls reproduce.
- Positions: `Slide-tags/coords/<lane>_coords.csv` (attached in step 3) and `andy_raw_data/<lane>/Positions/` (added
  in step 5).

To start from an existing Cell Ranger or CellBender run, export `CELLRANGER_RAW` or `CELLBENDER_H5` (see `config.sh`).

## Record

`REPRODUCE.md` describes the re-run of 6-8 Nov 2026. `record/` holds the scripts and logs of the run behind the
objects.
