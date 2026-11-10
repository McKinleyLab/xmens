# resources/annotations: how the cell-type labels were decided

Written 2026-11-09, after the fact, to state in one place what the July 2026 commits and the October 2026 public
export (`nhilgert/xmens-pipeline`) each record about the annotation. Nothing here changes a label.

## The record

- Decision tables, one row per Leiden-3.0 cluster (`cluster,label`): `resources/annotations/<species>_<lineage>.csv`;
  commits `7224f8e`, `2678268` (2026-07-29).
- Resolution, marker groups and label source per lineage: `resources/lineage_spec.yaml`; `7224f8e`, `0d184c9`,
  `2678268`.
- Coarse (whole-tissue) marker panels, Fig. 1c of the 2025 preprint: `resources/markers/panel_{mouse,human}.tsv`
  and its `README.md`; `0099f7e`, `bd7b036` (2026-07-21).
- The generator that produced the decision notebooks: `workflow/scripts/tools/gen_04c_notebooks.py`; first
  committed `8c84132` (2026-08-10), untracked in July on purpose.
- Review notebooks and rendered dot plots with fixed flag rules: `notebooks/review/*_review.ipynb`,
  `results/annotation_review/`; `73d702e` (2026-10-31).
- Per-cluster marker evidence (`note` column) and executed evidence notebooks: `nhilgert/xmens-pipeline`,
  `resources/annotations/*.csv` and `notebooks/annotation/`; `222bdc4`, `54dfa1c`, `95ee532`, `03c2714`
  (2026-10-26).

## How the labels were decided (2026-07-29)

For each lineage, the clustered object (`data/interim/latent_<species>_<lineage>.h5ad`, `leiden_scVI_3.0`) was
annotated in a notebook generated from `lineage_spec.yaml`. For every cluster the notebook showed one table with:

- the number of nuclei, doublet fraction and QC medians;
- `coarse_vote`: the most frequent 2025-submission lineage label among the cluster's nuclei, and its share;
- `fine_vote` (mouse stromal only): the most frequent 2025 subtype label (`cell_subtype_annotations` of
  `20251008_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad`, 97.1% of nuclei matched by barcode);
- the mean expression of each marker group in `lineage_spec.yaml` (for stromal: fibroblast, decidual,
  proliferation, contamination) and of every individual marker gene;

plus a marker dot plot by cluster and UMAPs of the lineage's `umap_genes`. A label was prefilled for each cluster
(auto-EXCLUDE on doublet evidence; otherwise the vote; for mouse epithelium a luminal/glandular marker rule) and
accepted or overridden from that table. The accepted labels are the committed CSVs. So the stromal subtypes
(eFib1, eFib2, eFib2i, eFib3, eFibP, dec1, dec2, dec3, decAp) are the 2025 subtypes carried onto the re-fit
clustering by per-cluster majority, checked against the lineage marker groups; they were not called from markers.
Agreement with the 2025 objects: mouse stromal 84% of nuclei keep their subtype, 99.99% their lineage
(`7224f8e`); human cell_type 92.80%, disagreements at the fibroblast/decidual boundary (4,052) and
stromal/smooth-muscle (413) (`2678268`). Mouse epithelium: luminal/glandular by markers (Wnt7a, Muc1, Lgr5 vs
Foxa2, Ihh, Lif). Immune fine labels: Alan Y. Baez Vazquez's per-nucleus table (`<species>_immune_fine.csv`).

The generated decision notebooks were not kept (`894e59b`, "regenerable"), so where a prefilled label was
overridden by hand is not recorded; the tables are the record of the outcome.

## The October 2026 evidence

On 2026-10-26 the public export regenerated each table with a `note` column giving the marker detection
fractions per cluster (`tools/write_annotation_notes.py`) and one executed evidence notebook per table
("No label is assigned in a notebook"), and removed the vote machinery and the 2025 label tables from that copy.
Its commit `222bdc4` describes the labels as "per-cluster marker calls"; the cluster-to-label pairs were unchanged
("byte-identical to the previous commit"). Read with this file: the notes are post-hoc marker evidence for labels
that were decided in July as described above. The October review (`notebooks/review/`, fixed flag rules) found
no change to make ("Changes against the v8 table of this clustering: none").

## What to cite

Methods: "Cluster labels were assigned per lineage at Leiden resolution 3.0. Stromal subtypes were carried over
from the 2025 analysis by per-cluster majority of the 2025 labels (84% of nuclei retained their subtype) and
checked against fibroblast, decidual and proliferation marker expression; epithelial compartments were assigned
from markers; immune subtypes from the per-nucleus annotation of A. Baez Vazquez. Decision tables and evidence:
`resources/annotations/`, `resources/lineage_spec.yaml`, `notebooks/review/`."
