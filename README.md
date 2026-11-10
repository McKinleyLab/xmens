# xmens: computational analysis for Cevrim et al.

*Induction of menstruation in mice reveals the regulation of menstrual shedding*, Cevrim et al.
(bioRxiv 2025.10.08.681007; revision in preparation, November 2026). Slide-tags single-nucleus spatial
transcriptomics of the mouse X-Mens model (six timepoints) and human endometrium (secretory and menstrual).

One folder per analysis step. Each folder holds the scripts that run the step, a `record/` of the run the paper
reports (inputs, checksums, logs, and `SOURCES.tsv` naming where each script came from), a comparison with the
previous run, and a README. Data objects are not in the repository; paths to them on the cluster are in each
folder's `config.sh` or README.

| folder | step | status (2026-11-09) |
|---|---|---|
| `preprocess/` | FASTQs -> Cell Ranger -> CellBender -> QC -> objects | record 2026-11-06; re-run identical |
| `annotation/` | cell-type and subtype labels (C. Cevrim) | objects 2026-11-07 (human), 2026-11-09 (mouse) |
| `pseudotime/` | harmonized time-informed pseudotime, fibroblast/decidual nuclei | in progress |
| `samap/` | cross-species mapping and pseudotime transfer | pending |
| `mi/` | pseudotime-associated genes (mutual information) | pending |
| `heatmaps/` | pseudotime heatmaps | pending |
| `de/` | differential expression (C. Cevrim) | pending |
| `scenic/` | regulon inference; removed from the paper 2026-11-09 | to be added (record 2026-11-08) |
| `notebooks/` | checks that read the objects and report, without writing | to be added |
| `docs/` | `HISTORY.md` (how the code got here), annotation notes | |
| `preprint_2025/` | the code of the October 2025 preprint (v1.0), unchanged | archive |

Earlier versions of this analysis: `preprint_2025/` (October 2025 notebooks) and the Snakemake pipeline of
June to November 2026 (`nhilgert/xmens_final`, to be imported here as `snakemake_2026/`). `docs/HISTORY.md`
explains the three and how to trace a script between them. The repository `nhilgert/xmens-pipeline` is an
archived simplified snapshot of the Snakemake version and should not be used.

Contributors: Nick Hilgert (pipeline, pseudotime, SAMap, MI), Cagri Cevrim (annotation, DE), Sean Eddy (code
review), Kara McKinley (PI).
