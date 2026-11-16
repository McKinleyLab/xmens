# environments

What each step ran with. The package lists were read from the environments on the FASRC cluster on 16 Nov 2026; no
package in these environments changed during the November runs. They record versions; they were not used to rebuild
the environments.

| environment or program | used by | record |
|---|---|---|
| `rsc_uv`: uv, Python 3.11.15 | `preprocess/` 3 (ingest) to 6, `pseudotime/`, `mi/`, `samap/` 0b-0c and 3-4, `pt_heatmaps/` 1 and 3, `notebooks/`, `figures/`, `tables/` | `rsc_uv.txt` |
| `samap_harmony_uv`: uv, Python 3.11.15 | `samap/` 1-2: standard SAMap (sc-samap 3.0.1); Harmony is installed but not used | `../samap/envs/packages_samap_harmony_uv.txt` |
| `cuda11.8_env`: conda, Python 3.10.13 | `preprocess/` 3 (QC, DoubletDetection 4.2) | `cuda11.8_env.yml` |
| `cellbender_022`: conda, Python 3.7.12 | `preprocess/` 2 (CellBender 0.2.2) | `cellbender_022.yml` |
| R 4.4.1 (FASRC module `R/4.4.1-fasrc01`) with Monocle 2.34.0 | `pt_heatmaps/` 2 | `R_monocle_sessionInfo.txt`; installed by `pt_heatmaps/0_install_monocle.sbatch` |
| Cell Ranger 7.2.0 (10x Genomics) | `preprocess/` 1 | |
| HMMER 3.4 (phmmer) | `samap/` 0a | |
| gcc 14.2.0 (FASRC module `gcc/14.2.0-fasrc01`) | `samap/`, for the libstdc++ that hnswlib needs | |
| `BulkRnaSeq`: conda | Çağrı Çevrim's DEG notebooks in `notebooks/` | `BulkRnaSeq.yml` |

`rsc_uv` and `samap_harmony_uv` use a Python that uv installed in the home directory of the account that built them
(`nhilgert`), so on the cluster they run only under that account. Elsewhere, create a Python 3.11.15 environment and
install the listed versions; the CUDA builds of torch need the PyTorch package index, and the RAPIDS packages
(`cudf-cu12` and others) NVIDIA's. `samap/envs/samap_env_pyproject.toml` is how `samap_harmony_uv` was specified;
its local `xmens` dependency was never installed, and the package list is the record.
