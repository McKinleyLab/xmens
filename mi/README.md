# mi: final objects -> MI with pseudotime -> human-mouse overlap -> supplementary table and Fig. 5F genes

Six steps, one script each, no Snakemake. Paths and settings live in `config.sh`. Submit from this folder; SLURM output
goes to `logs/`.

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/mi
export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_mi

sbatch 1_inputs.sbatch       # expression matrix and pseudotime per species   -> $OUT/inputs/{mouse,human}/
sbatch 2_mi.sbatch           # MI and 1,000 permutations, 1,700 array tasks   -> $OUT/chunks/{mouse,human}/
sbatch 3_combine.sbatch      # one table per species, with BH q               -> $OUT/mi_{mouse,human}.csv
sbatch 4_overlap.sbatch      # overlap of the top genes among orthologs       -> $OUT/overlap/
sbatch 5_table.sbatch        # the supplementary table (MI in bits)           -> $OUT/table/
sbatch 6_fig5f_genes.sbatch  # the 28 genes of the Fig. 5F heatmap            -> $OUT/fig5f/
```

Each step needs the one before it finished (steps 5 and 6 both need step 4). To queue all six at once:

```bash
j1=$(sbatch --parsable 1_inputs.sbatch)
j2=$(sbatch --parsable --dependency=afterok:$j1 2_mi.sbatch)
j3=$(sbatch --parsable --dependency=afterok:$j2 3_combine.sbatch)
j4=$(sbatch --parsable --dependency=afterok:$j3 4_overlap.sbatch)
sbatch --dependency=afterok:$j4 5_table.sbatch
sbatch --dependency=afterok:$j4 6_fig5f_genes.sbatch
```

- Step 2 is the long one: about 2,000 task-hours in the run of record (a mouse chunk takes 4-6 h, a human chunk about
  1 h), so roughly half a day with 200 tasks at a time. Tasks 0-199 are the mouse chunks, 200-1699 the human chunks.
- `bash 2_mi.sbatch --check` counts the finished chunks of each species without running anything (steps 1 and 3 also
  have `--check`).
- One chunk only: `sbatch --array=200 2_mi.sbatch` (human chunk 0), or `SLURM_ARRAY_TASK_ID=200 bash 2_mi.sbatch` in an
  interactive session.
- Outputs are never overwritten: a step or task whose output exists is skipped, and every output is written under a
  temporary name and renamed when complete. After a failure, submit the step again; only what is missing is made.
- Each Python script also runs on its own; `python scripts/mi_chunk.py --help` lists its arguments.

## Starting from earlier results

To skip steps 1-2 and start from the run of record's per-chunk results (the MI and p of every gene), export this before
submitting step 3:

```bash
export CHUNKS=/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/results/mi_dpt/pt_harmonized_{species}_nall/results
```

## Settings

- Nuclei: the fibroblast and decidual nuclei of the final objects with a pseudotime, `obs['pseudotime']` (rank
  normalized): 42,290 mouse (the 907 apoptotic decidual nuclei have none) and 47,533 human.
- Genes: detected in at least 1% of those nuclei and with at least 10 counts: 14,261 mouse, 15,008 human.
- Expression: log1p(CP10K) of the raw counts, each nucleus scaled over the tested genes.
- MI: scikit-learn 1.8.0 `mutual_info_regression`, 3 neighbours, `random_state` 0, on all nuclei; in nats in steps 2-4,
  in bits in the table (step 5).
- Permutations: 1,000 shuffles of the pseudotime per chunk, drawn by `numpy.random.default_rng(chunk)` and shared by the
  chunk's genes; p = (shuffles with MI >= observed + 1) / 1,001. 200 mouse and 1,500 human chunks (`config.sh`). The
  chunk count is part of the result, since each chunk draws its own shuffles.
- q: Benjamini-Hochberg over all tested genes of a species (statsmodels 0.14.6); significant: q < 0.05.
- Top genes: MI above the 90th percentile of the significant genes' MI (numpy, linear interpolation).
- Orthologs: Ensembl Compara 116 one-to-one protein-coding pairs (16,275). Overlap test: hypergeometric upper tail,
  `scipy.stats.hypergeom.logsf`; repeated at every percentile cutoff from the 40th to the 99th (Fig. S29C-D).
- Fig. 5F: pool = pairs in the top 10% by MI in both species among the pairs tested in both; per pool gene, Pearson r
  of the human and mouse z-scored profiles over 50 equal-count pseudotime bins; of the genes with r > 0.8, the 14 early
  (centroid < 0.4) and 14 late (> 0.6) genes with the highest r, in centroid order. The heatmap is drawn from these
  genes and the per-nucleus counts by the figure scripts, which are not part of this folder.

## Differences from the run of record

`record/` holds the unedited scripts, settings and small outputs of the run of record, with `SOURCES.tsv` (original paths
and sha256), `OUTPUTS.tsv` (the record's inputs and outputs, checksummed in place) and `LOGS_AND_OUTPUTS.md` (where that
run's logs are); `record/README.md` describes the run. Compared with it:

- Snakemake inputs became command-line arguments; the computations are unchanged.
- Step 1 reads the final objects and their `obs['pseudotime']`. The record read the 22 Oct 2026 objects and, for human,
  first copied the 30 Oct refit pseudotime into the object (`build_mi_record_input.py`). The nuclei, counts and
  pseudotime are the same, and so are the matrices (below).
- Both species run as one array. The record ran mouse as a 200-task array on 22-23 Oct and human as 1,500 Snakemake
  jobs on 31 Oct - 1 Nov.
- Step 3 stops if a chunk is missing; the record's merge warned and went on.
- Step 4 writes only what the Methods cite: the one-to-one pairs, the Methods overlap and the percentile sweep. The
  record's other overlap tables (two mouse control arms, a shuffled-pairing null, a 1-30% cutoff sweep and a
  shared-universe specificity table) are not ported.
- Step 5 is `build_human_mouse_mi_table.py` reorganized for reading. Its README sheet names step 5 and the input paths,
  and the label of the check against step 4 changed accordingly; every other cell is the same.
- Step 6 takes the MI tables and objects as files and no longer compares its selection with the drawn genes
  (`compare_with_record` does).
- Genes with equal q (most genes: p = 1/1,001) are listed in the order numpy's sort leaves them, which can depend on the
  CPU: on a Mac, the human table held the same rows in another order. No value depends on it; only the order of the
  four human genes with MI = 0 in the table's human_genes sheet follows it.

## Comparing with the run of record

```bash
sbatch compare_with_record.sbatch    # any time; compares the steps that have finished
```

The summary ends `logs/compare_with_record_<job ID>.out`; one row per comparison goes to
`$OUT/compare_with_record_<date>_<time>.tsv`. It compares step 1 with the record's matrices (`expr.npz`, `genes.npy`,
`pseudotime.npy`), each finished chunk with the record's chunk file, the combined tables with `mi_all.csv`, step 4 with
`overlap_of_record/`, the workbook with the approved workbook (`record/5_table/`) cell by cell, and step 6 with the
record's selection and the 28 genes drawn in Fig. 5F.

What reproduced in the checks of November 6, 2026 (`OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_xmens_mi_test`;
steps 3-6 started from the record's chunks with `CHUNKS`, and step 2 ran four chunks):

- Step 1: identical to the record's matrices: the same 14,261 mouse and 15,008 human genes, the same pseudotime of
  42,290 and 47,533 nuclei, and not one differing value in either expression matrix.
- Step 2: four chunks recomputed from the step-1 inputs, mouse chunks 0 and 199 (4.8 and 5.1 h) and human chunks 0
  and 1499 (57 and 61 min), byte-identical to the record's chunk files.
- Step 3: `mi_mouse.csv` and `mi_human.csv` byte-identical to the record's `mi_all.csv`.
- Step 4: the one-to-one pairs, `methods_overlap.json` (N 9,587, K 1,133, n 1,201, k 373, log10 p -83.7099) and the
  percentile sweep identical to `overlap_of_record/`.
- Step 5: every cell of the seven sheets identical to the approved workbook, apart from the build date, the input paths
  and one label (by design; the inputs' sha256 match).
- Step 6: the same pool of 273 pairs, the same 28 genes with the same r, and the same per-nucleus counts; the 28 genes
  are those drawn in Fig. 5F.
- On a Mac (numpy 2.4.6 and another math library), steps 3-5 gave the same values, but genes with equal q came out in
  another order, and log10 p differed in the 13th significant digit (-83.70989626027155 against -83.70989626027217).
