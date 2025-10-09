# xmens: cross-species single-cell analysis of induced menstruation

This repository accompanies the preprint *Induction of menstruation in mice reveals the regulation of menstrual shedding* by Çevrim *et al.*

---

**Version:** `v0.9` preprint release (October 2025).  
This snapshot includes figure generation notebooks and scripts for main and supplemental figures.  
A fuller release (`v1.0`) with processing and analysis pipelines will be available in the coming days.

---

### Figure index
- `fig3/fig_3a_spatial_human_cell_types.ipynb`: Fig. 3a  
- `fig3/fig_3c_pseudotime_heatmaps.ipynb`: Fig. 3c  
- `fig3/fig_3d_compare_pseudotime_gene_lists.ipynb`: Fig. 3d  
- `fig4/fig_4b_samap_menstrual.ipynb`: Fig. 4b (menstrual)  
- `fig4/fig_4b_samap_secretory.ipynb`: Fig. 4b (secretory)  
- `fig4/fig_4c_spatial_subtypes.ipynb`: Fig. 4c  
- `fig4/fig_4d_pseudotime_gradient.ipynb`: Fig. 4d  
- Selected supplemental notebooks under `supp/`

---

### Environment
Example conda environment for figure reproduction:  
`environments/sc_cuda_11.8.yml`

```bash
conda env create -f environments/sc_cuda_11.8.yml
conda activate sc_cuda_11.8

