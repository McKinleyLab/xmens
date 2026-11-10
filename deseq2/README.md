# deseq2: the paper's three pseudobulk DESeq2 fits, one notebook each

Each notebook rebuilds its pseudobulks from the object of record and fits PyDESeq2 with the settings of record. It then
checks its results against the table of record before it saves anything. All three reproduce their tables of record
exactly.

| notebook | fit | design | pseudobulks | genes tested | DE (up / down) |
|---|---|---|---|---|---|
| `xmens_decidual_vs_fibroblast.ipynb` | X-Mens decidual vs fibroblast nuclei, with the uninduced fibroblasts (27 Oct 2026) | `~ timepoint + group` | 11 | 16,784 | 4,703 (2,520 / 2,183) |
| `human_decidual_vs_fibroblast.ipynb` | our endometrium (Slide-tags, 2 donors), fit `B_new_celltype_release` (29 Oct) | `~ donor + group` | 4 | 20,211 | 3,198 (1,867 / 1,331) |
| `ulrich_decidual_vs_C7_fibroblast.ipynb` | Ulrich et al. 2024 atlas, decidual stroma vs C7-Fibroblast, contrast `dS_vs_C7fib_all_donors` (28 Oct) | `~ donor + state` | 20 (10 donors) | 13,772 | 886 (337 / 549) |

What the three share:
- **Test:** a Wald test of decidual vs fibroblast; a positive log2 fold change means higher in decidual cells.
- **Settings:** Cook's outlier refit on, independent filtering off.
- **Pseudobulks:** at least 20 nuclei (cells, for the atlas) each.
- **Genes tested:** at least 1 CPM in at least half of the decidual or of the fibroblast pseudobulks.
- **DE:** adjusted p < 0.05 and |unshrunk log2 fold change| > 1.

The human fit has 1 residual degree of freedom. Quote its fold changes and per-donor directions, not its p-values, as
evidence.

## Running

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/deseq2
export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_deseq2

sbatch xmens.sbatch     # 3-8 min on 8 CPUs, peak 3.0 GB   -> $OUT/xmens/
sbatch human.sbatch     # 2.5-6 min, peak 3.6 GB           -> $OUT/human/
sbatch ulrich.sbatch    # 1-2.5 min, peak 4.0 GB           -> $OUT/ulrich/
```

- **Order:** the three are independent, so they can be queued together.
- **Runtime:** each fit takes under 30 seconds. The rest is importing the packages and reading the object from lab
  storage, so times vary with file-system load. On 6 November, importing took 3 minutes in the two jobs that shared a
  node. The 30-minute Slurm limit leaves ample room.
- **Expected warning:** the human notebook prints a PyDESeq2 `UserWarning` that with fewer than 3 residual degrees of
  freedom the dispersion prior is poorly estimated. That is expected with 1 residual degree of freedom, and the
  results still match the record exactly.
- **Dry run:** `bash xmens.sbatch --check` (or either other script) lists what the run reads and writes, and whether
  each exists, without running anything.
- **Outputs:** a run is skipped when its output folder exists. The executed notebook goes to `$OUT/notebooks/`, with
  every number it printed, and the Slurm log goes to `logs/`.
- **Check failure:** a run that differs from the record stops at "Check against the record", before it saves anything.
- **By hand in Jupyter:** run `export JUPYTER_PATH=<this folder>/jupyter` before starting the server. Then open a
  notebook, choose the kernel `rsc` and run all cells. The default `OUT_DIR = None` writes nothing.

## Inputs

| notebook | object (read where it is) | copied into this folder |
|---|---|---|
| X-Mens | `xmens_final/data/interim/final_mouse.h5ad` | `resources/xmens_feature_manifest.tsv`; the table, pseudobulks and summed counts of record in `record/` |
| human | `xmens_final/data/deliverables/2026_09_15/Cevrim_Human_Only_Fib_Dec_Cells.h5ad` | the table of record in `record/` |
| Ulrich | `/n/netscratch/eddy_lab/Everyone/nhilgert/ulrich2024/ulrich2024_20donor_stromal.h5ad` | `resources/ulrich_feature_manifest.tsv`; the table and pseudobulks of record in `record/` |

- **Copied files:** `record/` and `resources/` hold unchanged copies from `xmens_final/resources/ivd/`, taken on
  6 November 2026. `SHA256SUMS` lists them.
- **The Ulrich atlas is on netscratch,** which deletes files older than 90 days (60 from the December maintenance). The
  copy is dated 18 September 2026. It is an unmodified CELLxGENE download (dataset `c93622aa…`), and the notebook
  checks its SHA-256. To get it back:
  1. Run the notebook by hand once with `H5AD` set to a lab-storage path and `DOWNLOAD_H5AD = True`. A fresh download
     was checked on 6 November: same SHA-256, same results.
  2. Point `ULRICH_H5AD` in `config.sh` at that path.
- **Release objects:** the fits of record were made on the objects above, not on the release objects in
  `final_objects/`.
  - The human release object gives identical pseudobulks. Run with `-p H5AD <path>`: it has no donor column, so the
    notebook takes the donor from the sample names. The check then differs only in the `object` column, which names
    the input file.
  - The X-Mens release objects hold the same counts, but their labels and timepoints are renamed (`cell_type`,
    `1.5_dpi_early`, …), so the settings cell would need editing.

## Verification (6 November 2026)

- **Installed folder:** the three sbatch scripts were run with `sbatch` from this folder (Slurm jobs 58391127-58391129,
  run folder `/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_deseq2`). All three completed and stated
  "identical to the record".
- **Tables:** each results table is byte-identical to its table of record: the same genes, every column and every DE
  call.
- **Pseudobulks:** identical to the record's (nuclei and summed counts).
- **Environment:** `rsc_uv` (Python 3.11.15, PyDESeq2 0.5.4, NumPy 2.2.0).
  - The fits of record ran in `rsc_uv_py31111` (Python 3.11.11, the same PyDESeq2 and NumPy). The `rsc_uv` runs here
    reproduce them byte for byte, and the Ulrich notebook also gives identical output in `rsc_uv_py31111`.
  - The X-Mens and Ulrich fits are also identical on Intel Cascade Lake, Intel Ice Lake and AMD Genoa nodes.
  - With other NumPy versions the human p-values move: NumPy 2.4.6 gives 3,195 DE genes instead of 3,198. Keep this
    environment.
- **Differences from the code of record** (`DE_IVD_code_2026-10-31`), none of which changes a result:
  - **X-Mens:** the code of record sums features that share a gene. None do, so the notebook asserts that instead. It
    prints the inputs' SHA-256 rather than requiring them, so other objects can be run through it.
  - **Human:** a fallback, used only when the object has no donor column, that takes the donor from the sample names.
  - **Human:** only the primary table is reproduced. The sensitivity and regression fits are left out: `coarse_refined`
    labels, one pseudobulk per Slide-tags section, and the old notebook's settings. The section fit is the lane check
    in the results note of 29 October.
  - **Ulrich:** only the reported contrast is run. The release's other contrasts and its IVD summaries are left out.

## Files

- `config.sh`: paths and environment, sourced by every sbatch script. It keeps all caches out of the home folder.
- `xmens.sbatch`, `human.sbatch`, `ulrich.sbatch`: one per notebook.
- `jupyter/kernels/rsc/kernel.json`: the `rsc_uv` kernel, the same one as in `samap`.
- `resources/`, `record/` and `SHA256SUMS`: see Inputs.
- `logs/`: the Slurm logs.
