# scenic: fibroblast/decidual objects -> 100 SCENIC runs per species -> consensus regulons -> Fig. 5E, Fig. S28

Six steps, one script each, no Snakemake. Paths and settings live in `config.sh`; the 200 runs (seeds 1000-1099 in
each species) are the rows of `runs.tsv` (row n = array task n of steps 2 and 3; rows 1-100 mouse, 101-200 human).
Submit from this folder; SLURM output goes to `logs/`.

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/scenic
export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_scenic

sbatch 1_looms.sbatch              # pooled objects, 1% gene filter        -> $OUT/scenic_seeds/<species>/<species>_dsc_pooled.*
sbatch 2_grnboost2.sbatch          # 200 runs: 15k sketch + GRNBoost2      -> .../looms_sketch15k/, adjacencies_sketch15k/
sbatch 3_cistarget.sbatch          # 200 runs: cisTarget pruning           -> .../regulons_sketch15k/
sbatch 4_consensus_aucell.sbatch   # consensus regulons + AUCell           -> .../consensus_sketch15k_weighted_soft_occ6*/
sbatch 5_pseudotime_rho.sbatch     # rho with pseudotime, Fig. 5E table    -> $OUT/correlations/
sbatch 6_figures.sbatch            # Fig. 5E and Fig. S28                  -> $OUT/figures/
```

Each step needs the previous one finished. To queue all six at once:

```bash
j1=$(sbatch --parsable 1_looms.sbatch)
j2=$(sbatch --parsable --dependency=afterok:$j1 2_grnboost2.sbatch)
j3=$(sbatch --parsable --dependency=afterok:$j2 3_cistarget.sbatch)
j4=$(sbatch --parsable --dependency=afterok:$j3 4_consensus_aucell.sbatch)
j5=$(sbatch --parsable --dependency=afterok:$j4 5_pseudotime_rho.sbatch)
sbatch --dependency=afterok:$j5 6_figures.sbatch
```

Step 2 is most of the time: GRNBoost2 took 4 to 7.5 hours per run on 16 cores in September 2026. Steps 3-6 are
minutes per task.

- `bash 2_grnboost2.sbatch --check` (every step has it) prints each task's inputs and whether its outputs exist,
  without running anything.
- One run only: `sbatch --array=101 2_grnboost2.sbatch` (row 101 of runs.tsv: human, seed 1000), or
  `SLURM_ARRAY_TASK_ID=101 bash 2_grnboost2.sbatch` in an interactive session. Steps 1 and 4 have one task per
  species (1 mouse, 2 human).
- Outputs are never overwritten: a task whose output exists is skipped, and files and folders are written under a
  temporary name and renamed when complete.
- Each Python script also runs on its own; `python scripts/ctx_only.py --help` lists its arguments.

## What each run records

Every Python call goes through `scripts/run_with_settings.py`, which runs the script unchanged and writes
`$OUT/logs/<log name>.settings.json` next to its log (for example `3_cistarget_mouse_seed1000.settings.json`):

- `options`: the script's options as its own parser resolved them, defaults included (`options_from_defaults` names
  those not given on the command line); `script_constants`: its module-level constants (for example the cisTarget
  databases of `ctx_only.py`).
- `calls`: every argument of the library functions that do the computation, with the library's defaults filled in
  (`library_defaults` names those not passed): `gs_gap` (the sketch), `modules_from_adjacencies`, `prune2df` and
  `df2regulons` (cisTarget), `aucell` (the consensus script's AUCell and `pyscenic aucell`).
- `subprocesses`: for GRNBoost2, arboreto's options, the call it makes for every target gene, and every parameter of
  the gradient-boosting regressor (scikit-learn), early stopping included.
- `inputs` (size, date, sha256 up to 256 MB), the Python version and every installed package's version, thread and
  SLURM variables, start and end times, and the exit code.

The file is written before the computation starts and updated before each call, so a run that fails or is killed
still shows what it ran with. Until the call succeeds it is named `<name>.settings.json.partial.<pid>`. `pyscenic
aucell` parses `--rank_threshold` and `--nes_threshold` but does not use them; `calls` shows what `aucell` received.

## Starting from the run of record

The 100 runs per species of September 2026 can stand in for steps 2 and 3. Set these before submitting:

```bash
# its sketches and GRNBoost2 adjacencies (then run steps 1 and 3-6)
export GRN_FROM=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/scenic_v2_2026_08_28
# its cisTarget regulons as well (then run steps 1, 4, 5 and 6)
export REGULONS_FROM=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/scenic_v2_2026_08_28
```

The shortest full check of the deterministic steps, about an hour:

```bash
export REGULONS_FROM=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/scenic_v2_2026_08_28
j1=$(sbatch --parsable 1_looms.sbatch)
j4=$(sbatch --parsable --dependency=afterok:$j1 4_consensus_aucell.sbatch)
j5=$(sbatch --parsable --dependency=afterok:$j4 5_pseudotime_rho.sbatch)
j6=$(sbatch --parsable --dependency=afterok:$j5 6_figures.sbatch)
sbatch --dependency=afterok:$j6 compare_with_record.sbatch
```

## Settings

- pySCENIC 0.12.1 (arboreto 0.1.6, ctxcore 0.2.0) in `.conda/envs/pyscenic` for steps 1-4; `rsc_uv` (scanpy 1.11)
  for steps 5-6. Both are named in `config.sh`.
- Step 1 input: the fibroblast/decidual objects of release `2026_09_15`, raw counts from `layers['counts']`. The
  pooled loom keeps genes detected in at least 1% of nuclei (14,245 of 33,698 mouse, 15,011 of 36,601 human); it is
  the expression matrix of GRNBoost2, cisTarget and the AUCell of record.
- Step 2: each run trains GRNBoost2 (pySCENIC's `arboreto_with_multiprocessing.py`, 16 workers) on its own
  15,000-nucleus geometric sketch of the pooled loom's genes, drawn on `obsm['X_scVI']`; the run's seed seeds both the
  sketch and GRNBoost2. TF lists `allTFs_mm.txt` and `allTFs_hg38.txt`.
- Step 3: pySCENIC's Python API with its defaults: co-expression modules from positive correlations without dropout
  masking, `prune2df` with rank threshold 1,500, AUC threshold 0.05 and NES 3; motif collection v10 (clustered),
  mm10 or hg38 rankings at 500 bp up/100 bp down and 10 kb up/down of the TSS.
- Step 4: a TF is kept if it formed a regulon in at least 80 of the 100 runs and keeps at least 10 targets; a target
  is kept if its link to the TF recurs in at least 6 runs, weighted by the fraction of the TF's runs that contain it.
  This gives the regulons of record: 193 mouse, 141 human. The same at a floor of 2 runs feeds Fig. S28. Activity of
  record: `pyscenic aucell --weights` on the pooled loom, ranking seed 0, AUC threshold 0.05.
- Step 5: Spearman rho of each regulon's activity with `obs['pseudotime']` of the final objects (the 907 apoptotic
  mouse decidual nuclei have none and are left out: 42,290 mouse and 47,533 human nuclei). TFs are paired through
  Ensembl 116 one-to-one orthologs: 111 TFs have a regulon in both species.
- Step 6: Fig. 5E colours a pair when |rho| >= 0.2 in both species and names 13 TFs; Fig. S28 shows the TF and target
  thresholds.

## Differences from the run of record

`record/` holds the scripts that made the run of record (unedited copies), `SOURCES.tsv` (original paths and sha256
of its inputs, scripts and outputs), `LOGS_AND_OUTPUTS.md` (where its logs are), `step_counts.tsv` (counts at every
step) and `outputs/` (the record's Fig. 5E table, rho tables and figures, which `compare_with_record` checks against).
Compared with them:

- The Snakemake rules and the one-off jobs of October became six scripts with command-line arguments; the
  computations are unchanged.
- Every Python call runs through `scripts/run_with_settings.py`, which writes the settings files described above.
  The run of record wrote none; its settings are in its logs, manifests and the scripts in `record/scripts/`.
- Step 1 reads release `2026_09_15`. The run of record read the netscratch `2026_08_28_*` objects on 8 Sep 2026,
  which were later rewritten with added `obs` columns; the release is a copy of the rewritten files. Nuclei, genes,
  raw counts, `timepoint` and `X_scVI` are the same in both.
- `build_consensus_aucell.py` takes the lane folder as a required `--rev` (the record's copy defaulted to a retired
  path); nothing else in the script changed.
- Step 4 builds, in one job, what the record built three times: floor 2 with the Snakemake rule (9-10 Sep), floor 6
  with a one-off job (14 Oct) and the AUCell of record (26 Oct). The YAML export uses the script the repository has
  since adopted (`scenic_regulons_to_yaml.py`) instead of the inline code of `cliauc_run.sh`; both call pySCENIC's
  `save_to_yaml` and check the reload.
- Step 5 runs the 2 Nov 2026 version of the correlation script, which reads `obs['pseudotime']` of the final objects.
  The table of record was made on 30 Oct by a one-off script; the 2 Nov version reproduced it byte for byte.
- Fig. S28: the manuscript file was drawn in a Claude Science session (`record/scripts/regulon_recurrence_cutoffs_fig.py`,
  Helvetica Neue). Step 6 draws it with the repository's port of that script, which writes the same text and numbers
  in the cluster's fonts. The floor-2 link tables it reads were made in September by a script that was not kept;
  step 6 writes them with `s28_consensus_edges.py` from the floor-2 consensus.

## Comparing with the run of record

```bash
sbatch compare_with_record.sbatch    # any time; compares the steps that have finished
```

The summary ends `logs/compare_with_record_<job ID>.out`; one row per comparison goes to
`$OUT/compare_with_record_<date>_<time>.tsv`. It compares the pooled objects (genes, nuclei, counts, scVI latent),
each run's sketch and adjacencies, each run's regulons, both consensus floors and the AUCell of record with the
record's files, the rho tables and the Fig. 5E table with `record/outputs/`, the Fig. 5E SVG byte for byte with the
manuscript panel, and Fig. S28's text with the manuscript file.

What reproduced in the checks of November 6, 2026 (run folder `2026_11_06_xmens_scenic`; steps 1 and 4-6 run from
the record's cisTarget regulons, steps 2-3 for seed 1000 of each species):

- Step 1: the pooled objects are identical to the record's (genes, nuclei, counts, scVI latent).
- Steps 2 and 3, seed 1000 of each species: the sketches contain the same nuclei and genes, GRNBoost2 gave the same
  links and importances (1.23 M and 1.24 M links; 5.6 and 6.2 h on 16 cores), and cisTarget the same regulons as the
  record's.
- Step 4: both consensus floors and the AUCell of record are byte-identical to the record's files.
- Step 5: rho of every regulon is identical, and the Fig. 5E table is byte-identical (18 TF pairs up in both species,
  8 down in both, 14 up in human and down in mouse, none down in human and up in mouse).
- Step 6: the Fig. 5E SVG is byte-identical to the manuscript panel. Fig. S28 has the same 84 text labels and numbers
  as the manuscript file, and its two input tables are identical.
- `run_with_settings.py` changes no result. cisTarget, `pyscenic aucell`, a sketch, GRNBoost2 (on a 300-gene test
  loom), the floor-6 consensus and the correlations gave the same output through it as without it, and as the record
  (`tmp/2026_11_06_scenic_settings_test/`).
