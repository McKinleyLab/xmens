# environments

| environment | used by | packages |
|---|---|---|
| `rsc_uv` (uv, Python 3.11.15) | most steps | `rsc_uv.txt` |
| `samap_harmony_uv` (uv, Python 3.11.15) | `samap/` steps 1-2 (SAMap 3.0.1; Harmony not used) | `../samap/envs/` |
| `cuda11.8_env` (conda) | `preprocess/` step 3 (QC) | `cuda11.8_env.yml` |
| `cellbender_022` (conda) | `preprocess/` step 2 (CellBender 0.2.2) | `cellbender_022.yml` |
| R 4.4.1, Monocle 2.34.0 | `figures/pt_heatmaps/` | `R_monocle_sessionInfo.txt` |
| `BulkRnaSeq` (conda) | `Mouse DEGs` and the DEG analyses `01`-`04` in `notebooks/` | `BulkRnaSeq.yml` |

Also Cell Ranger 7.2.0 (`preprocess/` step 1), HMMER 3.4 (`samap/` step 0a) and gcc 14.2.0 (`samap/`).
