# mi: mutual information between expression and pseudotime, and the human-mouse overlap

For each species, the mutual information (MI) of every gene with pseudotime and a permutation p-value; then the
overlap of the top genes among one-to-one orthologs, the supplementary table and the genes of Fig. 5F. One script
per step. Paths and settings are in `config.sh`.

## Run

```bash
cd mi && mkdir -p logs
export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_mi    # a new folder for each run

j1=$(sbatch --parsable 1_inputs.sbatch)                              # expression and pseudotime per species
j2=$(sbatch --parsable --dependency=afterok:$j1 2_mi.sbatch)         # MI and 1,000 permutations, 1,700 array tasks
j3=$(sbatch --parsable --dependency=afterok:$j2 3_combine.sbatch)    # one table per species, with BH q
j4=$(sbatch --parsable --dependency=afterok:$j3 4_overlap.sbatch)    # overlap of the top genes among orthologs
sbatch --dependency=afterok:$j4 5_table.sbatch                       # the supplementary table (MI in bits)
sbatch --dependency=afterok:$j4 6_fig5f_genes.sbatch                 # the genes of the Fig. 5F heatmap
```

Step 2 is the long one: a mouse chunk takes 4-6 h, a human chunk about 1 h. `bash 2_mi.sbatch --check` counts the
finished chunks. A step or task whose output exists is skipped, so after a failure submit the step again.

## Settings

- Nuclei: the fibroblast and decidual nuclei with a pseudotime. Genes: detected in at least 1% of them, with at least
  10 counts.
- Expression: log1p(CP10K) of the raw counts, each nucleus scaled over the tested genes.
- MI: scikit-learn 1.8.0 `mutual_info_regression`, 3 neighbors, `random_state` 0.
- Permutations: 1,000 shuffles of pseudotime per chunk (200 mouse and 1,500 human chunks); p = (shuffles with MI >=
  observed + 1) / 1,001. q: Benjamini-Hochberg over the tested genes of a species; significant at q < 0.05.
- Top genes: MI above the 90th percentile of the significant genes' MI.
- Orthologs: Ensembl Compara 116 one-to-one protein-coding pairs, matched to the objects by Ensembl gene ID. Overlap
  test: hypergeometric upper tail.
- Fig. 5F: among pairs in the top 10% by MI in both species, Pearson r of the human and mouse z-scored profiles over
  50 equal-count pseudotime bins; of the genes with r > 0.8, the 14 early (centroid < 0.4) and 14 late (> 0.6) genes
  with the highest r.

## Record

`record/RUNS.md` lists the runs behind the paper. `record/` also holds the scripts and logs of earlier runs.
