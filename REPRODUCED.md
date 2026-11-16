# Reproduction

Each step below was run a second time from code in this repository, and its outputs were compared with those of the
run behind the paper (the run of record). "Committed as" means the rerun used the working tree and was committed
right after; the committed files are the ones that ran. Run logs, job IDs, checksums and comparison tables are
archived in commit `3605739` (`<step>/record/`); the comparisons of 15-16 Nov are in
`/n/lab_storage/eddy_lab/users/nhilgert/xmens_reproduced/` on the FASRC cluster.

| step | run of record | rerun | result |
|---|---|---|---|
| Preprocessing, `preprocess/` | Cell Ranger Feb 2026; CellBender and QC Jul 2026 | 14-15 Nov, clean clone at `09e3d27` | identical at every step |
| Pseudotime, mouse | 10 Nov, `pseudotime_mouse_final_1110.ipynb` at `e5a12fc` | 12 Nov, the same notebook | object identical (h5diff) |
| Pseudotime, human | 12 Nov, `pseudotime_human_final_1112.ipynb` at `0b83615` | 12 Nov, the same notebook | object identical (h5diff) |
| MI steps 1-3 | 12 Nov at `33a8f5c` | 12 Nov at `33a8f5c`, both species from scratch | MI tables identical |
| MI steps 4-6 | 13 Nov at `4bbaf54` and `0f407dd` | 15 Nov, clean clone at `v2.0` | no result differs |
| MI table, values only | 13 Nov, step 5 at `0f407dd` | 15 Nov at `3051518` | same values; no formulas |
| SAMap step 1 | 12 Nov at `33a8f5c` | 12 Nov at `33a8f5c` | identical |
| SAMap steps 0a-0c and 2-4 | 14 Nov at `19929ae` | 14 Nov at `19929ae` | homology tables and every result table identical |
| Heatmaps, `pt_heatmaps/` | 14 Nov | 15 Nov, committed as `404c480` | every figure and z matrix identical |
| Figures, `figures/` | 13-14 Nov | 15 Nov, committed as `16f2dc3` | every figure file identical or equal in content |
| QC table, `tables/` | 14 Nov | 16 Nov, committed as `3605739` | data sheets identical |

Notes:

- Code since the reruns: each step's code in this commit is that of its rerun. For SAMap, comments and docstrings
  changed after `19929ae`; no code line did.
- Pseudotime: the reruns were copies of the notebooks with only the output path changed. h5diff found the objects
  identical except `uns/pseudotime/computed`.
- MI steps 1-3: the run of record reused the mouse chunks of 11 Nov; the rerun computed both species again. Of 1,728
  files, 1,715 were identical, 5 equal in content (time stamps and paths) and 7 were logs. The one other file was that
  run's workbook, whose README sheet differed; steps 4-6 were redone on 13 Nov.
- SAMap: the fit of 12 Nov did not reproduce and was replaced by the runs of 13 and 14 Nov. In the rerun of 14 Nov,
  the 19 result tables of steps 0c-5 were byte-identical (step 5, the figures, has since moved to `figures/`), and
  `Samap_Final.h5ad` was equal in 41 of its 42 datasets; the other, `uns/provenance`, holds time stamps and the run's
  path.
- Heatmaps: of 35 files, 28 identical, 2 equal in content (a gz file after decompression; the Monocle PDF without its
  dates), and the logs and executed notebook not compared. `inputs/cells_human.tsv.gz` differs: it carries the final
  human subtype names, which the heatmaps do not use.
- Figures: of 89 figure files, 55 identical and 34 equal in content (SVG without dates and ids, PDF without creation
  date). Fig. S32 was then redrawn with new labels (`f1094eb`); only panels d, f and g changed.
- QC table: rebuilt from the deposited objects with a shorter README sheet.

Not rerun: the annotation (step 2), the notebooks that added columns to the deposited objects (step 6, run at
`d7db734` on 14 Nov), and the renaming of `eFib` to `Fib` in the human deposited objects (16 Nov; compared with the
objects before it, `cell_type` is the only field that changed). Çağrı Çevrim's other notebooks are not covered here.

## To rerun

Run the steps in the order of the table in `README.md`, each into a new output folder (`OUT`); each step's README
gives the commands, and the notebooks stop rather than overwrite an output.

- Environments: `rsc_uv` and `samap_harmony_uv` run only under the account that built them; see
  `environments/README.md`. The others can be used by members of the eddy_lab group.
- Inputs on scratch (`/n/netscratch`, purged after 90 days): `figures/S30_mi_gene_selection.ipynb` and
  `pt_heatmaps/config.sh` read the `mi/` runs of 15 and 13 Nov; `tables/` reads the QC run of 6 Nov; the pseudotime
  notebooks read a reference object for a comparison at their end. Rerun the upstream step and point them at its
  output, or use the copies made on 16 Nov in `/n/lab_storage/eddy_lab/users/nhilgert/xmens_scratch_copies/`.
- The QC table also reads `samples.tsv` and the QC `config.yaml` from
  `/n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/preprocess/`. The first equals `preprocess/samples.tsv`; the
  second is `preprocess/record/config/config.yaml` in commit `3605739`.
- Preprocessing step 3 runs on AMD nodes (`--constraint="genoa|turin"`), on which the doublet calls reproduce. The
  pseudotime notebooks fit scVI on a GPU.
