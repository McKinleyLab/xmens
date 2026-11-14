# samap: cross-species alignment of the human and mouse fibroblast/decidual nuclei

SAMap on the final fibroblast/decidual objects, then mouse cell type and pseudotime transferred onto the human
nuclei, statistics and figures. One notebook per step in `notebooks/`, run by its sbatch script with papermill.
Paths are in `config.sh`.

## Run

```bash
cd samap && mkdir -p logs
unset PHMMER_HITS MAPS SOURCE_MAPS SAMAP_INPUTS SAMAP_JOINT SAMAP_FINAL
export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_samap    # a new folder for each run

a=$(sbatch --parsable --array=1-11,24-41,43-55,58-64 0a_phmmer.sbatch)       # phmmer, each proteome against the other
t=$(sbatch --parsable --array=12-23,42,56,57 --mem=64G 0a_phmmer.sbatch)     # the chunks with titin need more memory
b=$(sbatch --parsable --dependency=afterok:$a:$t 0b_homology.sbatch)         # best hit per gene pair
i=$(sbatch --parsable 1_inputs.sbatch)                                       # raw counts of the final objects
export SOURCE_MAPS=$OUT/0_homology
c=$(sbatch --parsable --dependency=afterok:$b:$i 0c_relabel.sbatch)          # gene names as in the objects
export MAPS=$OUT/0_homology_relabeled
s=$(sbatch --parsable --dependency=afterok:$c 2_samap.sbatch)                # the SAMap fit, about 1 h
x=$(sbatch --parsable --dependency=afterok:$s 3_transfer.sbatch)             # cell type and pseudotime, mouse to human
y=$(sbatch --parsable --dependency=afterok:$x 4_statistics.sbatch)           # tables, edge shares, permutation control
sbatch --dependency=afterok:$y 5_figures.sbatch                              # Fig. 6C, 6E and the SAMap supplement
```

`bash <step>.sbatch --check` prints what a step reads and writes, without running anything. A step whose output
exists in `OUT` is skipped.

## Outputs

| step | writes in `$OUT` |
|---|---|
| 0a-0c | `0_homology/phmmer/`, `0_homology/hsmm/`, `0_homology_relabeled/hsmm/` |
| 1 | `1_inputs/{hs,mm}_fib_dec_all_raw.h5ad` |
| 2 | `2_samap/samap_joint.h5ad`, `samap.pkl` |
| 3 | `3_transfer/Samap_Final.h5ad`, `transfers.csv` |
| 4 | `4_statistics/` |
| 5 | `5_figures/` |

The executed notebooks, with everything they printed, are in `$OUT/notebooks/`.

## Environments

Kernels in `jupyter/kernels/`: `samap` (sc-samap 3.0.1, anndata 0.12.19, scanpy 1.11.5, hnswlib 0.8.0) for steps 1
and 2, `rsc` (anndata 0.12.11, scanpy 1.11.5) for the others. phmmer is HMMER 3.4. Package lists are in
`record/config/`.

## Record

`record/RUNS.md` lists the runs behind the paper. `record/` also holds the scripts and logs of earlier runs.
