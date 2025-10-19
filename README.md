This repository accompanies the preprint   *Induction of menstruation in mice reveals the regulation of menstrual shedding* by Çevrim *et al.* (2025).

bioRxiv: [https://doi.org/10.1101/2025.10.08.681007](https://doi.org/10.1101/2025.10.08.681007)

---

**Version:** `v1.0` preprint release (October 2025).  
This snapshot includes scripts and notebooks to reproduce main figures and analyses.
Pipeline pieces are organized by module (processing, SAMap, pseudotime).

---

### Repository structure

```text
xmens/
├── environments/           # conda environments
├── processing/             # QC, integration, clustering
├── samap/                  # cross-species mapping
├── pseudotime/             # pseudotime
├── fig3/                   # Figure 3 notebooks
├── fig4/                   # Figure 4 notebooks
├── supp/                   # supplemental notebooks
└── README.md
```

---

### Environments
Minimal conda environment for figure reproduction:  
`environments/sc_cuda_11.8.yml`

```bash
conda env create -f environments/sc_cuda_11.8.yml
conda activate sc_cuda_11.8
```

For running cross-species SAMap integration scripts, use the `samap_env.yml` environment.
