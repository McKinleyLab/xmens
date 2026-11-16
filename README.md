# xmens

Code for Çevrim *et al.*, *Induction of menstruation in mice reveals the regulation of menstrual shedding*
([bioRxiv](https://doi.org/10.1101/2025.10.08.681007)).

| | folder | |
|---|---|---|
| 1 | `preprocess/` | FASTQs to QC'd objects |
| 2 | `notebooks/Integration, Cleaning, and Annotation.ipynb` | annotation |
| 3 | `pseudotime/` | pseudotime |
| 4 | `mi/` | genes associated with pseudotime |
| 5 | `samap/` | human-mouse alignment (SAMap) |
| 6 | `notebooks/` | the final objects, with pseudotime, the SAMap transfers and the Stro-7 score |
| 7 | `figures/pt_heatmaps/` | heatmaps along pseudotime |
| 8 | `figures/` | figures |
| 9 | `tables/` | QC table |

`notebooks/` also holds the Stro-7 scoring, the mouse DEG comparison with oil-induced decidualization and the
CC3/pH3 segmentation. `orthology/` holds the ortholog table used in step 4; `environments/` lists the software.

Data will be available upon publication.
