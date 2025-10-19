# Pseudotime analysis

This directory contains notebooks, scripts, and SLURM job files used to compute diffusion pseudotime and to estimate mutual information (MI) between gene expression and pseudotime.

---
## Directory structure

```text
pseudotime/
├── notebooks
│   ├── integated_human_pseudotime.ipynb
│   ├── integrated_mouse_pseudotime.ipynb
│   └── unintegrated_mouse_pseudotime.ipynb
└── scripts
    ├── mi_sbatch.sh
    └── mi_worker.py
```

---

## Overview

### Pseudotime inference notebooks

The notebooks in `notebooks/` perform diffusion pseudotime (DPT) analysis for the human and mouse datasets (integrated and unintegrated versions). Unintegrated DPT operates on all mouse fibroblast cells, while integrated DPT is computed on a subset of human and mouse cells at comparable timepoints along the (induced) menstrual cycle. MI analysis is performed on integrated DPT.

### Mutual information (MI) analysis

Scripts in `scripts/` compute mutual information between gene expression and pseudotime across expressed genes.
- `mi_worker.py` computes MI for one gene chunk and permutation-based _p_-values  
- `mi_sbatch.sh` SLURM array launcher for parallel MI computation

**Example usage**

Run the MI pipeline on a SLURM-based compute cluster:
```text
cd pseudotime/scripts
sbatch --array=0-999 mi_sbatch.sh
```
Each job runs one chunk (of 1000 chunks) using `mi_worker.py` and saves output to `results/mi_###.csv`.

