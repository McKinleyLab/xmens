# samap: proteomes -> homology maps; final fib/dec objects -> SAMap -> transfers -> statistics and figures

Steps 0-5, no Snakemake. Step 0 makes the homology maps SAMap needs: phmmer searches (0a, an array of shell tasks),
then a notebook (0b). Every later step is one notebook in `notebooks/`, run by its sbatch script with papermill; the
executed copy, with every number the step printed, is the step's log. Paths live in `config.sh`. Submit from this
folder; SLURM output goes to `logs/`, the executed notebooks to `$OUT/notebooks/`.

To check the science, read the notebooks' first code cell after the parameters: in step 2 it is the SAMap fit
(one SAM per species, `SAMAP(...)`, `sm.run(...)`, as in the SAMap tutorial; 49 lines), in step 3 both transfers
(76 lines, half of them loading the graph and the labels).
Everything after a "Save" or "Bookkeeping" heading only writes files, evaluates, or records checksums. The sbatch
scripts are plumbing: find `config.sh`, check inputs, skip if done, run the notebook.

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/samap
mkdir -p logs                # SLURM writes the jobs' output here and does not create the folder
export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_samap_full  # a new folder for each run

sbatch 0a_phmmer.sbatch      # phmmer, each proteome against the other     -> $OUT/0_homology/phmmer/
sbatch 0b_homology.sbatch    # homology maps from the phmmer hits          -> $OUT/0_homology/hsmm/
sbatch 1_inputs.sbatch       # raw counts of the final objects             -> $OUT/1_inputs/
sbatch 2_samap.sbatch        # the SAMap fit (about 1 h, 8 cores, 128 GB)   -> $OUT/2_samap/
sbatch 3_transfer.sbatch     # cell type and pseudotime, mouse to human     -> $OUT/3_transfer/Samap_Final.h5ad
sbatch 4_statistics.sbatch   # tables, edge shares, permutation control     -> $OUT/4_statistics/
sbatch 5_figures.sbatch      # Fig. 6C, 6E, cell identity, Supp. Fig. S30   -> $OUT/5_figures/
```

Each step needs the one before it finished; step 1 does not need step 0. To queue everything at once:

```bash
j0a=$(sbatch --parsable 0a_phmmer.sbatch)
j0b=$(sbatch --parsable --dependency=afterok:$j0a 0b_homology.sbatch)
j1=$(sbatch --parsable 1_inputs.sbatch)
j2=$(sbatch --parsable --dependency=afterok:$j0b:$j1 2_samap.sbatch)
j3=$(sbatch --parsable --dependency=afterok:$j2 3_transfer.sbatch)
j4=$(sbatch --parsable --dependency=afterok:$j3 4_statistics.sbatch)
sbatch --dependency=afterok:$j4 5_figures.sbatch
```

- `bash 2_samap.sbatch --check` (any step) prints what the step reads and writes and whether each exists, without
  running anything.
- Outputs are never overwritten: a step whose output exists is skipped. Steps 1-4 write into `$OUT/<step>.partial/`
  and rename it when complete; step 5 skips a figure whose executed notebook exists and redraws one that did not
  finish. A step stops instead of skipping when its existing output was made from other inputs (a start-from
  variable changed between runs into one `OUT`), so use a new `OUT` for each run.
- Before reading them, each step checks the files from outside this folder (proteomes, final objects) against the
  sha256 checksums in `config.sh`, and stops if one differs.
- The notebooks of steps 0b-4 run by hand too: in Jupyter, with `export JUPYTER_PATH=<this folder>/jupyter` before
  the server starts, open one, choose kernel `samap` (steps 1-2) or `rsc` (steps 0b, 3, 4), set `OUT` in the
  parameters cell and run all cells; the other parameters default to the step before's output in `OUT`, as in
  `config.sh` (the sbatch scripts also check input checksums, which a by-hand run does not). The samap kernel needs
  `module load gcc/14.2.0-fasrc01` in the shell that starts Jupyter. The figure notebooks of step 5 take their tables
  and output folder as parameters; run them through `5_figures.sbatch`.

| step | notebook | kernel | reads | writes |
|---|---|---|---|---|
| 0a | (shell: `awk` and `phmmer`, 64 array tasks) | | the two proteomes | `0_homology/chunks/`, `0_homology/phmmer/{hu_vs_mu,mu_vs_hu}/*.tbl` |
| 0b | `0_homology.ipynb` | rsc | step 0a, the two proteomes | `0_homology/hsmm/{hs_to_mm,mm_to_hs}.txt`, `homology.json` |
| 1 | `1_inputs.ipynb` | samap | the final fib/dec objects | `1_inputs/{hs,mm}_fib_dec_all_raw.h5ad`, `inputs.json` |
| 2 | `2_samap.ipynb` | samap | step 1, the homology maps | `2_samap/samap_joint.h5ad`, `samap.pkl`, `eval.json`, `mapping_*.csv`, `figures/`, `meta.json` |
| 3 | `3_transfer.ipynb` | rsc | step 2, the final objects | `3_transfer/Samap_Final.h5ad`, `transfers.csv`, `summary.json` |
| 4 | `4_statistics.ipynb` | rsc | step 3, the final objects | `4_statistics/tables/`, `graph_block/`, `edge_share/`, `alignment_scores.csv`, `label_agreement.csv`, `ptnull/`, `statistics.json` |
| 5 | `5_figures/*.ipynb` (11) | rsc | step 4 | `5_figures/fig6/`, `5_figures/sfig_samap/` (panels in `panels/`) |

## Starting from earlier results

Export one of these before submitting. The later steps read it instead of the earlier steps' output, and the steps
whose output would go unused skip themselves, so all steps can still be submitted (or queued as above):

```bash
# the phmmer hit tables of record (step 0a skips)
export PHMMER_HITS=/n/lab_storage/eddy_lab/users/nhilgert/data/proteomes/hm
# the homology maps of record (steps 0a and 0b skip)
export MAPS=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/samap_integrate/homology/hsmm_hmmer
# the inputs of record (step 1 skips)
export SAMAP_INPUTS=/n/lab_storage/eddy_lab/users/nhilgert/code/samap-scvi/full_fibdec/inputs_raw
# the fit of record, as exported on 2 Sep 2026 (steps 0-2 skip, including the 1 h fit)
export SAMAP_JOINT=/n/lab_storage/eddy_lab/users/nhilgert/code/samap-scvi/full_fibdec/runs_raw/vanilla_raw/samap_transferred.h5ad
# the deposited object (steps 0-3 skip: statistics and figures only)
export SAMAP_FINAL=/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects/Samap_Final.h5ad
```

## Settings

- Homology maps: phmmer (HMMER 3.4) of every human protein (Ensembl release 113 GRCh38 pep.all, 123,845 proteins)
  against every mouse protein (release 113 GRCm39 pep.all, 66,558) and back, hits kept at E <= 1e-5, in 32
  round-robin chunks per direction, 4 threads each. For each (query gene, target gene) pair, the best bit score over
  their isoforms and that hit's E-value; genes are the `gene_symbol` of the FASTA headers.
- Inputs: the raw counts (layer `counts`) of the deposited fib/dec objects in `final_objects/`: 47,533 human nuclei
  (36,601 genes) and 43,197 mouse nuclei (33,698 genes).
- SAM, for each species (SAMap's defaults, sc-samap 3.0.1): counts scaled to the median library size, log2(x + 1),
  values of 1 or less set to 0, genes expressed in more than 96% of nuclei left out; 3,000 genes, 100 principal
  components, 20 neighbors.
- SAMap: three iterations with the hmmer human-mouse homology maps (`hsmm_hmmer`), 8 cores. Seeded since 12 Nov 2026
  (`SEED = 0` in `notebooks/2_samap.ipynb`, with hnswlib held to one thread); the fit of record predates this (see
  "What reproduced").
- Transfers, over the human x mouse block of the joint graph: cell type by edge-weighted vote of the mouse
  neighbors; pseudotime as the edge-weighted mean of the mouse pseudotime of record (rank / N), over mouse neighbors
  that have one (decAp nuclei have none).
- Statistics: edge shares with the final objects' cell types, the mean total cross-species edge weight per nucleus,
  and SAMap's alignment scores between compartments (computed from the graph as SAMap computes them); agreement of
  the transferred cell type with the human labels; permutation control with 1,000 shuffles per design (seed
  20261022), Spearman rho against the human pseudotime of record and Moran's I over the 10 nearest nuclei within each
  section.
- Figures: the notebooks of record with two edits (see "Differences from the run of record"): fixed 0-1 viridis,
  values as stored, every nucleus drawn.

## The Methods, claim by claim

Where each statement of the Methods section "Cross-species alignment and transfer of categorical and continuous
values (SAMap)" is made or measured. The numbers match a run started from the fit of record (`SAMAP_JOINT`); a fresh
fit gives slightly different ones (see "Comparing with the run of record").

| Methods | here |
|---|---|
| Ensembl release 113 proteomes; phmmer (HMMER 3.4) in both directions, E <= 1e-5 | `config.sh`; `0a_phmmer.sbatch`. Release checked: `record/methods_check/ensembl_release_check.txt` |
| hits collapsed to gene symbols (`gene_symbol` field), best isoform pair per ordered symbol pair; proteins without a symbol excluded | step 0b; counts in `0_homology/hsmm/homology.json` |
| SAMap keeps pairs with E <= 1e-6 (its default) | inside SAMap (`samap/core/homology.py`): it keeps pairs found in both directions at E <= 1e-6, weights each by the mean bit score of the two directions, and drops pairs below a quarter of the best pair of both their genes |
| sc-samap 3.0.1 on raw integer counts indexed by gene symbol; SAMap's internal settings; three iterations; 20 cross-species neighbors; edge weights rescaled by expression correlation | step 1 (inputs); step 2: `PREPROCESS` and `SAM_RUN` equal the settings SAMap uses when it builds the SAMs itself; `n_iterations=3`; `cross_species_k=20` and `scale_edges_by_corr=True` are `SAMAP.run` defaults; `2_samap/meta.json` |
| edge shares 88%, 82%, 85%, 92% | step 4, `edge_share/samap_compartment_edge_share.csv` |
| mean total cross-species edge weight 12.6-13.5 of at most 20 | step 4, `edge_share/samap_compartment_edge_weight.csv` |
| alignment scores 0.58 and 0.58 matched, 0.08 and 0.12 off-diagonal | step 4, `alignment_scores.csv` (with the final labels; step 2's `mapping_celltype.csv` gives the same for a fresh fit) |
| the transfer operator; cell type and pseudotime from mouse to human; 47,509 and 47,352 human nuclei with a value | step 3, `transfers.csv`, `summary.json` |
| transferred vs human labels; Spearman rho vs human pseudotime, all and per compartment | step 4, `label_agreement.csv`; the observed row of `ptnull/ptnull_statistics_table.csv` |
| permutation control (42,290 mouse nuclei, two designs, 1,000 each, Spearman and Moran's I) | step 4, `ptnull/`; Supplementary Fig. S30G in step 5 |
| human to mouse: Leiden clusters and Stro-7 cosine; 43,196 of 43,197 mouse nuclei with a value | not in this pipeline: xmens_final `workflow/rules/41_fig_c7_presence.smk` (`c7_tables.py`), from the same joint object. The one mouse nucleus without a human edge is counted in step 4's `samap_compartment_edge_weight.csv` |

## Environments

The two kernels in `jupyter/kernels/` start the Python of an existing environment; `config.sh` points papermill at
them with `JUPYTER_PATH`.

- `samap`: `samap_harmony_uv` (sc-samap 3.0.1, anndata 0.12.19, scanpy 1.11.5, hnswlib 0.8.0); needs
  `module load gcc/14.2.0-fasrc01` (`config.sh` loads it).
- `rsc`: `rsc_uv` (papermill, anndata 0.12.11, scanpy 1.11.5, matplotlib 3.10.8).
- phmmer: HMMER 3.4 built in `/n/lab_storage/eddy_lab/users/nhilgert/software/hmmer-3.4` (`config.sh`, `PHMMER`).

`record/config/` has the package lists of both as installed on 6 Nov 2026, which describe the environments that ran
(rebuild from these), and, for reference, the spec and lock file of xmens_final's `workflow/envs/samap`, which pins
numpy < 2 and is not the environment step 2 runs in.

## Differences from the run of record

`record/` holds the unedited scripts, job files and logs that made the SAMap results of the paper, with
`SOURCES.csv` (original paths and sha256), `OUTPUTS.csv` (the large inputs and outputs, checksummed in place),
`numbers.csv` (every number the figures and methods cite) and `LOGS_AND_OUTPUTS.md`. Compared with them:

- One notebook per step replaces the scripts, sbatch files and Snakemake rules; the computations are unchanged.
- Step 0a makes its own chunks, with one `awk` line (round robin); the record's chunks came from a splitter nobody
  recorded, and are byte-identical to these. It adds `--noali` to phmmer, which keeps the alignments out of the
  summary output (1 GB per chunk in the record) and leaves the hit tables unchanged. The phmmer of record was the
  one on the PATH of 24 Apr; this one is HMMER 3.4 as rebuilt on 4 Aug 2026.
- Step 1 reads the final objects; the record read the netscratch `2026_08_28` objects. Their counts, nuclei and
  genes, in order, are identical (checked 6 Nov 2026). The inputs carry the final objects' names and labels in obs
  and no embeddings; the fit reads neither.
- Step 2 writes the joint object itself, before the evaluation; the record exported it from `samap.pkl` a day later,
  with the label transfers attached. The evaluation's label-based numbers use the final objects' labels.
- Step 3 votes with the final mouse `cell_type` and transfers the final mouse `pseudotime`, joined by exact nucleus
  name. The record voted with the `2026_08_28` labels and transferred `pt_harmonized` of `2026_10_22_flagged_harmonized`,
  rank-normalized; the 2 Nov checks found the same transferred cell types, and the final `pseudotime` is that
  rank-normalized value. It writes only the two transfers `Samap_Final.h5ad` keeps (the record's transfer object also
  had the raw and within-sample pseudotime variants and two more label transfers) and writes `Samap_Final.h5ad`
  directly (the record made the transfer object on 24 Oct and cut its obs on 2 Nov).
- Step 4 reads the final objects and `Samap_Final.h5ad` directly; the record went through plot tables of the same
  objects and checked the human pseudotime against `human_pt_2026_10_30_run4_a40/cells.parquet` (equal on every
  nucleus, 3 Nov 2026). `human_plot_table.parquet` holds the human object's obs only (no embeddings or gene columns).
- The `embedding_outlier` flags are gone (decided 6 Nov 2026). No panel has hidden the flagged nuclei since 1 Nov
  2026, so the figures do not change; the record's tables and figure provenance records carry the flags (6 human and
  46 mouse nuclei), this pipeline's do not. Their tables sit in `final_objects/old/embedding_outliers/`.
- Step 5 runs the figure notebooks of record (xmens_final `59e00dd`, `notebooks/samap/`) with two changes: they no
  longer read or count the `embedding_outlier` flags, and the cell-identity notebook imports `cross_species_panels.py`
  from beside it instead of from the repository. Their opening cells still describe running them in xmens_final.

## Comparing with the run of record

```bash
sbatch compare_with_record.sbatch    # any time; compares the steps that have finished
```

The summary ends `logs/compare_with_record_<job ID>.out`; one row per comparison goes to
`$OUT/compare_with_record_<date>_<time>.csv`. It compares the proteome chunks, phmmer hit tables and homology maps
with those of 24 Apr 2026, the inputs with `inputs_raw/`, the fit with
`runs_raw/vanilla_raw` (alignment scores, graph edges, UMAP), `Samap_Final.h5ad` with the deposited one (obs column
by column, every other dataset by checksum), the step-4 tables and statistics value by value and the figures pixel by
pixel with the 3 Nov 2026 set (`record/4_5_figures/results_final/`; SVGs apart from their drawing date and element
ids), and the numbers of `record/numbers.csv` that the run computes. Inputs taken from the run of record are "not
compared"; step 0 is "incomplete" until all 32 chunks and hit tables of a direction exist.

What reproduced in the checks of 6 and 7 November 2026 (run folders on netscratch, each with its
`compare_with_record_<date>_<time>.csv`):

- Homology maps (step 0b, from the phmmer hit tables of record; `PHMMER_HITS` as above): both maps byte-identical to
  the maps of record (12,833,563 and 12,780,006 hits; 1,615,662 and 1,616,943 gene pairs; the five known orthologs
  present). Step 0a, tested on chunk 1 of each direction: the chunks are byte-identical to the record's and the hits
  the same in every field; in one table the target-name column is padded one character narrower for 18 queries,
  which the maps do not see. The full 64-task search has not been repeated (about 40 min on 64 x 4 cores).
- The SAMap runs below read the maps of record (`MAPS` was fixed to them until step 0 was added).
- Inputs (step 1): counts, nuclei and genes identical, in order, to the inputs of record.
- Starting from the fit of record (`SAMAP_JOINT` as above; `2026_11_06_xmens_samap_from_fit_of_record`), steps 3-5
  reproduce the paper. `Samap_Final.h5ad` equals the deposited object in every obs column, and in every dataset
  outside obs except `uns/provenance`. The one difference: the transferred pseudotime differs from the stored values
  by at most 1.1e-16, a rounding difference in the last bit (the permutation control of record found the same when it
  recomputed the transfer). The transferred cell type is identical on all 47,533 human nuclei. The step-4 tables,
  graph block, edge shares and permutation statistics are identical to the 3 Nov set (the transferred-pseudotime
  column of `joint.parquet` carries the same 1.1e-16 difference), and so are the 66 numbers of `record/numbers.csv`
  the run recomputes. All 23 figure PNGs are pixel-identical, the 15 SVGs differ only in their
  drawing date and element ids, and the one CSV is identical.
- A fresh fit (`2026_11_06_xmens_samap`) does not give the same graph: on identical inputs and settings, 72% of the
  record's graph edges are in the new graph (Jaccard index 0.57 over all edges; for a human nucleus's mouse
  neighbors 0.43 on average, the same set for 229 nuclei), the UMAP layout differs, and the mean alignment score is
  0.650 (record 0.660). The results the paper reports move little:
  - transferred pseudotime: Pearson r 0.989 and Spearman rho 0.987 with the record; a value on 47,392 human nuclei
    (record 47,352)
  - transferred cell type: the same as the record's on 46,370 of 47,533 human nuclei; agreeing with the human labels
    on 83.4% of labelled nuclei (record 84.3%)
  - Spearman rho with the human pseudotime: 0.735, 0.722 and 0.273 for all, decidual and fibroblast nuclei (record
    0.734, 0.728, 0.260); Moran's I 0.66, 0.549 within decidual and 0.339 within fibroblast (record 0.66, 0.553,
    0.338)
  - permutation control: every observed statistic is above all 1,000 shuffles of both designs, except fibroblast rho
    against within-compartment shuffles, which 4 of 1,000 reach (record 6)
  - edge shares: within 0.018 of the record's
  - figures: drawn on a different UMAP layout, so they differ from the paper's pixel by pixel
- After the review fixes of 7 Nov 2026 (`2026_11_07_xmens_samap_*`): from the fit of record, with every step queued
  (steps 0-2 skipped themselves), steps 3-5 reproduce the paper as above, and the 87 numbers the run computes equal
  `record/numbers.csv`, including the edge-weight totals and alignment scores of the Methods. From the deposited object
  (`SAMAP_FINAL`), step 4 and 5 outputs are identical. A second fresh fit, on maps that step 0b made from the hit
  tables of record (byte-identical to the maps of record), again gave its own graph: 88% of the record's edges
  (Jaccard 0.78; mouse neighbors of a human nucleus 0.62), mean alignment 0.665. Transferred pseudotime: Pearson r
  0.992 with the record, a value on 47,384 human nuclei; transferred cell type the same on 46,739; agreement with the
  human labels 83.9%; rho 0.732, 0.723, 0.263; Moran's I 0.641, 0.542, 0.321; permutation control as in the record
  (fibroblast rho against within-compartment shuffles reached by 7 of 1,000); edge shares within 0.015; alignment
  scores 0.58 and 0.59 matched, 0.08 and 0.13 off-diagonal; mean total cross-species edge weight 12.65-13.67.

To reproduce the paper's numbers and figures exactly, start from the fit of record (`SAMAP_JOINT`); a full run from
the final objects reproduces them to the precision above. The run-to-run differences come from SAMap drawing random
numbers from numpy's unseeded global state (SAMap issue #135) and from hnswlib building its neighbor indexes on
several threads: since 12 Nov 2026 the fit cell seeds numpy and holds every hnswlib index to one thread, and two fits
on the same inputs are then identical bit for bit (joint graph, UMAP, alignment scores; test in
`$NS/2026_11_12_samap_seed_test/`, two seeded fits against each other and against an unseeded re-run). The fit of
record of 10 Nov 2026 was unseeded and cannot be regenerated exactly.

Run times on 6 Nov 2026: step 0a 40 min per task (64 tasks; 50-55 min in the run of record), step 0b 20 min,
step 1 2 min, step 2 56 min, step 3 1-2 min, step 4 2 min, step 5 4-6 min, `compare_with_record` 1-3 min.
