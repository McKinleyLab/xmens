# In vivo vs in vitro decidualization: differential expression comparisons

Notebooks comparing gene expression changes of decidualization in vivo (single-cell / single-nucleus
RNA-seq, pseudobulked) with in vitro decidualization models (bulk RNA-seq), in human and mouse.

## Notebooks

| # | Notebook | What it does |
|---|---|---|
| 01 | `01_DEG_Analysis_Human_InVitro.ipynb` | Bulk RNA-seq of human endometrial stromal fibroblasts decidualized with cAMP + MPA (3, 6, 8 days): FASTQ -> FastQC/MultiQC -> kallisto -> PyDESeq2; defines the in vitro core DEGs (same direction at all three timepoints) |
| 02 | `02_DEG_Analysis_Human_InVivo_vs_InVitro.ipynb` | Pseudobulk DE of decidual vs non-decidual stromal cells in three human single-cell datasets (this study, Ulrich et al. 2024, HECA), compared with the in vitro core DEGs |
| 03 | `03_DEG_Analysis_Mouse_InVivo_vs_InVitro.ipynb` | Mouse: pseudobulk DE of decidual vs non-decidual stromal cells in X-Mens; in vitro decidualization DEGs (Wang et al. 2020) compared with natural pregnancy (Wang et al. 2020), artificial decidualization (Wang et al. 2020) and X-Mens |
| 04 | `04_DEG_Analysis_Mouse_InVivo_Model_Comparisons.ipynb` | Mouse: pseudobulk DE of X-Mens decidual vs non-decidual stromal cells, then natural pregnancy (Wang et al. 2020), artificial decidualization (Wang et al. 2020) and X-Mens DEGs each classified in the other two |

Run order: 02 reads `all_genes_merged.tsv` written by 01. Notebooks 03 and 04 are independent and can be
run in any order.

Each notebook has the same layout: Libraries, Datasets, Variables, analysis steps, plots.
Paths are placeholders (`/directoryto/...`); set them in the *Datasets* / *Variables* sections.
After changing a setting (e.g. the animals or excluded subtypes), run all cells again.

## Methods summary

* **Quantification (bulk):** kallisto 0.52.0, single-end mode (fragment length 200 +/- 20), against the
  Ensembl 110 human transcriptome (cDNA + ncRNA), matching Cell Ranger `GRCh38-2024-A`; transcript
  counts and TPM summed to genes.
* **Pseudobulk (single-cell):** raw counts summed per biological sample (donor / animal) x group;
  pseudobulks with < 20 cells dropped. Paired design (`~donor + condition`) when samples contain both
  groups, otherwise `~condition`.
* **Biological samples vs 10x lanes:** in this study's data, `timepoint` identifies the biological sample
  and `sample` the 10x lanes. Lanes are summed into their sample and never treated as replicates
  (human: 2 donors x 3 lanes; mouse: one animal per timepoint, 1-2 lanes each). In notebooks 03 and 04
  the animals used and any excluded cell subtypes are set in section 3.
* **Differential expression:** PyDESeq2 0.5.4 (Wald test, Benjamini-Hochberg adjusted p-values;
  apeglm-shrunken log2FC reported alongside).
* **DEG cutoffs (all analyses):** padj <= 0.01, |fold change| >= 2 (|log2FC| >= 1), and mean
  expression > 3 TPM (bulk) or CPM (pseudobulk) in at least one group - matching Wang et al. 2020
  (fold change >= 2, FDR <= 0.01).
* **Gene matching:** Ensembl gene IDs (version suffix removed); for Wang et al. 2020 (mm9 annotation),
  genes whose ID is not found are matched by gene symbol. The Wang et al. tables list DEGs only, so a gene
  absent from a table counts as "not a DEG" in that model.
* **Comparisons:** each reference list is shown as a stacked bar (shared up / shared down / shared,
  opposite direction / reference DEG only). Notebook 04 also reports the same-direction agreement among
  shared DEGs and the agreement expected between unrelated lists,
  p_up(A) * p_up(B) + p_down(A) * p_down(B).

### Note on GSE114296 (Rytkonen et al. 2019)

Three SRA runs are linked to the wrong GEO sample. Decidual marker expression (IGFBP1, FOXO1, PRL) and
the authors' own count table (GEO supplementary file) show that runs labelled ESF_1, ESF_2 and
DSC3days_3 are decidualized (3 days) and ESF_3, DSC3days_1 and DSC3days_2 are undifferentiated controls.
Notebook 01 uses the corrected assignment.

## Data

Raw and processed data are not included in this repository.

| Dataset | Source |
|---|---|
| Stadtmauer & Wagner 2022, bulk RNA-seq (6 d cAMP + MPA vs base media) | SRA PRJNA679324: SRR13081985, SRR13081986, SRR13081989, SRR13081990 |
| Rytkonen et al. 2019, bulk RNA-seq (3 d / 8 d decidualization) | GEO GSE114296 (SRA SRP145267): SRR7142614-SRR7142620, SRR7142622 |
| Ulrich et al. 2024, stromal cells, 20 donors | `Ulrich_OnlyStroma_20Donors.h5ad` (authors' public release) |
| Mareckova et al. 2024, Human Endometrial Cell Atlas (HECA) | `endometriumAtlasV2_cells_with_counts.h5ad` (authors' public release) |
| Wang et al. 2020, mouse DEG tables | Supplementary Tables S2 (NPD), S3 (AD), S6 (IVD) - place in `BulkRnaSeq/Wang2020/` |
| This study (human and mouse snRNA-seq) | `*_final.h5ad` |

For the Rytkonen and Stadtmauer FASTQs, download the runs above (e.g. from ENA as `.fastq.gz`) into
`BulkRnaSeq/Rytkonen2019/` and `BulkRnaSeq/Stadtmauer2022/`. The ENA files for SRR7142620 and
SRR7142622 already contain the second sequencing runs (SRR7142621, SRR7142623). Notebook 01 downloads
the Ensembl 110 transcriptome and builds the kallisto index if they are not present.

## References

* Stadtmauer DJ, Wagner GP (2022). Single-cell analysis of prostaglandin E2-induced human decidual cell
  in vitro differentiation: a minimal ancestral deciduogenic signal. *Biol Reprod* 106(1):155-172.
  https://doi.org/10.1093/biolre/ioab183
* Rytkonen KT et al. (2019). Decidualization of human endometrial stromal fibroblasts is a multiphasic
  process involving distinct transcriptional programs. GEO GSE114296.
* Ulrich ND et al. (2024). Cellular heterogeneity and dynamics of the human uterus in healthy
  premenopausal women. *PNAS*. https://doi.org/10.1073/pnas.2404775121
* Mareckova M et al. (2024). An integrated single-cell reference atlas of the human endometrium.
  *Nat Genet*. https://doi.org/10.1038/s41588-024-01873-w
* Wang et al. (2020). Comparative analysis of mouse decidualization models at the molecular level.
  *Genes* 11(8):935. https://doi.org/10.3390/genes11080935
* Muzellec B et al. (2023). PyDESeq2: a python package for bulk RNA-seq differential expression analysis.
  *Bioinformatics* 39(9):btad547. https://doi.org/10.1093/bioinformatics/btad547

## Environment

```bash
mamba env create -f environment.yml
mamba activate BulkRnaSeq
python -m ipykernel install --user --name BulkRnaSeq
```

Then select the `BulkRnaSeq` kernel in Jupyter.
