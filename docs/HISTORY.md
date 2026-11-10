# How the analysis code got here

Three generations of code have analyzed the X-Mens / human endometrium Slide-tags data. This page says what each
one is, why the current one replaced the previous one, and how to trace a script from one to the next. Written
2026-11-09, the day the group decided on the current layout; the per-step lines at the end are filled in as each
step is finished this week.

## 1. October 2025: the preprint code (`preprint_2025/`)

Notebooks and scripts that made the figures of the bioRxiv preprint (2025.10.08.681007): `processing/qc.ipynb`
and `integrate_and_cluster.ipynb`, SAMap, pseudotime, the Fig. 3, Fig. 4 and supplementary notebooks. Released
as `nhilgert/xmens` v1.0 (two commits, 2025-10-09 and 2025-10-19). Moved unchanged into `preprint_2025/` on
2026-11-09; its git history is attached to the files.

## 2. June to November 2026: the Snakemake pipeline (`snakemake_2026/`)

`xmens_final`, a Snakemake workflow (63 commits, 2026-06-29 to 2026-11-06) that re-ran the analysis from the
CellBender counts: QC, scVI integration, a CSV-driven annotation stage, pseudotime, SAMap, mutual-information
genes, SCENIC, zone modules and the figure panels. Its QC reproduced the preprint's nuclei to within 0.1%
(commit `1197108`); its annotation (2026-07-29, commits `7224f8e`, `0d184c9`, `2678268`) carried the 2025
labels onto a re-fit clustering by per-cluster majority, checked against marker-group expression
(`resources/annotations/README.md` there). It produced the release objects of 2026-09-15 through 2026-10-30 and
the figure drafts in circulation until 2026-11-09. A simplified public snapshot, `nhilgert/xmens-pipeline`
(14 commits, 2026-10-25/26), is archived; one of its commit messages describes the annotation as marker calls,
which overstates it (see the README above). The full repository is imported here as `snakemake_2026/` with
its history and is not run anymore.

## 3. November 2026: one script per step (this repository, top level)

On 2026-11-06 the preprocessing was rewritten as plain SLURM scripts (`preprocess/`), and on 2026-11-09 the group
(K. McKinley, C. Cevrim, S. Eddy, N. Hilgert) decided that the paper's analysis would be re-done this way from
the objects C. Cevrim re-annotated on 2026-11-06/09: one folder per step, each with a reviewed script, a kept run
of record, a comparison with the previous run, and a paragraph of documentation; everything pushed to this
repository as it is finished. SCENIC was removed from the paper at the same meeting.

Why the change. Not because the Snakemake version was wrong: its QC and its scripts are what the new folders
run. The group needed each step to be reviewable by a second person as a single file with its parameters
visible, and the annotation to be re-done rather than carried over. A DAG with config indirection and five
rounds of stage renumbering was hard to review that way, and the Snakemake annotation was anchored to the 2025
labels. One script per step, with the inputs frozen, was the smallest structure that met both needs.

What is reproducible across the change. The preprocessing re-run from FASTQs (2026-11-08) reproduced the
Snakemake pipeline's QC'd objects exactly (`preprocess/REPRODUCE.md`): Cell Ranger, CellBender, QC and doublet
calls identical in 15 of 15 lanes, objects byte-identical. The annotated objects used from here on are
`20261109_mouse_allcells_final.h5ad` and `20261107_Human_allcells_final.h5ad` (C. Cevrim), checked against the
QC'd objects in `notebooks/final_objects_check.ipynb`.

## How to trace a script

Every folder at the top level that took code from the Snakemake pipeline has a `record/SOURCES.tsv`: for each
file, its original path in `xmens_final`, the commit, the modification time and the sha256 of the copy. The
file is written when the folder is made and checked with `cmp` against the originals. `preprocess/record/`
is the model (34 files from commit `59e00dd`, all checked 2026-11-08). The line for each step below says
whether its code is unchanged, wrapped (same script behind a new sbatch front with the parameters at the top),
or rewritten.

| step | folder | from Snakemake | state | run of record | compared with |
|---|---|---|---|---|---|
| preprocessing | `preprocess/` | ingest, qc, qc_report, concat; rules 03-06 | wrapped | 11-06, rerun 11-08 | same |
| annotation | `annotation/` | none (new, C. Cevrim) | new | 2026-11-09 | 2025 labels and July tables (to fill) |
| pseudotime | `pseudotime/` | `pseudotime/harmonized_pseudotime.py`, rule 17b | to fill | to fill | Sep 15 record |
| SAMap | `samap/` | `scripts/samap/*`, rules 19, 21 | to fill | to fill | Oct 2026 run |
| MI genes | `mi/` | `scripts/mi/mi_pseudotime.py`, rule 22 | to fill | to fill | Sep 15 run of record |
| heatmaps | `heatmaps/` | `scripts/figures/...` (to fill) | to fill | to fill | Oct 2026 panels |
| DE | `de/` | `scripts/de/*`, rule 18 | to fill (C. Cevrim) | to fill | to fill |
| SCENIC | `scenic/` | `scripts/scenic/*`, rule 24 | wrapped | 2026-11-08 (mouse) | Sep run: core regulons hold |
