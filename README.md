# xmens

Code for Çevrim *et al.*, *Induction of menstruation in mice reveals the regulation of menstrual shedding*
([bioRxiv 2025.10.08.681007](https://doi.org/10.1101/2025.10.08.681007)), as revised in November 2026. Tag `v1.0` is
the code of the preprint.

Single-nucleus RNA-seq with spatial positions (Slide-tags) of mouse uteri through induced menstruation (X-Mens) and
of human endometrium, from the FASTQs to the figures. Each step is a folder with a README, its paths in `config.sh`
and one Slurm script or notebook per stage. The paths point to the FASRC cluster at Harvard, where the analysis ran.

## Steps

| | where | reads | writes |
|---|---|---|---|
| 1 | `preprocess/` | FASTQs | QC'd objects, one per species |
| 2 | `notebooks/Integration, Cleaning, and Annotation.ipynb` | step 1 | annotated objects |
| 3 | `pseudotime/` | step 2 | the fibroblast/decidual objects with pseudotime |
| 4 | `mi/` | step 3 | MI of each gene with pseudotime, the human-mouse overlap, the MI table, the Fig. 5F genes |
| 5 | `samap/` | step 3 | the SAMap object; mouse cell type and pseudotime transferred to the human nuclei |
| 6 | `notebooks/*_20261114.ipynb` | steps 2, 3 and 5 | the deposited fibroblast/decidual objects |
| 7 | `pt_heatmaps/` | steps 4 and 6 | Fig. 4H, Fig. 5F and the MI heatmaps (Monocle 2) |
| 8 | `figures/` | steps 4-6 | one notebook per figure |
| 9 | `tables/` | steps 1, 2 and 6 | the QC table |

`orthology/` holds the human-mouse ortholog table (Ensembl 116) used in step 4.

## Software

`environments/README.md` lists every environment and program, its version and the steps that use it. Most steps run
in `rsc_uv` (Python 3.11.15); SAMap steps 1 and 2 run in `samap_harmony_uv` (standard SAMap, sc-samap 3.0.1).

## Data

Raw reads and objects: [accession to be added]. The deposited objects are the all-cells objects of step 2 and the
fibroblast/decidual objects of step 6, with one later change: in the human objects the cell type `eFib` was renamed
`Fib` (Ç. Çevrim, 16 Nov 2026), and nothing else changed. On the cluster the steps read the objects from
`/n/lab_storage/eddy_lab/Lab/mckinley/data/`, which, like the environments, needs membership in the eddy_lab group.

## Reproduction

Steps 1, 3-5 and 7-9 were run a second time from committed code and their outputs compared with those behind the
paper. `REPRODUCED.md` lists the reruns, what the comparisons found, what was not rerun and how to rerun. The run
logs, job IDs, checksums and comparison tables are in commit `3605739`, in the folders named `record/`.

## Other notebooks (Çağrı Çevrim)

In `notebooks/`: `Integration, Cleaning, and Annotation.ipynb` (step 2), `Figures_Panels.ipynb`,
`Human-Stro7_scoring.ipynb`, `Mouse DEGs.ipynb` (X-Mens decidual vs fibroblast DEGs compared with oil-induced
decidualization), `CC3_pH3_Segmentation.ipynb` (CC3+ and pH3+ cells in the tissue images) and the bulk RNA-seq DEG
analyses `01`-`04_DEG_Analysis_*.ipynb`. `environments/BulkRnaSeq.yml` is the environment of the DEG notebooks.

## Citation

Çevrim Ç. *et al.* Induction of menstruation in mice reveals the regulation of menstrual shedding. bioRxiv (2025).
https://doi.org/10.1101/2025.10.08.681007
