# xmens_scenic: the SCENIC run of record (copy of record, 2026-11-06)

Unedited copies of the scripts that produced the SCENIC results in the manuscript: the regulons and activities in
`/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/scenic_v2_2026_08_28/scenic_seeds/<species>/`
(the "lane" below), the Fig. 5E table, Fig. 5E and Fig. S28. `SOURCES.tsv` gives each file's original path, size and
sha256 (each copy here was checked against its original); `step_counts.tsv` the counts at every step, read from the
record's own files; `LOGS_AND_OUTPUTS.md` where the logs are. Not runnable from here: paths inside the scripts point to
their original locations. The runnable version is the folder above this one.

| Step | Files | Ran | Output |
|---|---|---|---|
| 1 Pooled objects | scripts/scenic_make_loom.py; rule scenic_pooled (config/13_scenic.smk) | 2026-09-08 | lane `<species>_dsc_pooled.{h5ad,loom,manifest.json}` |
| 2 Sketch + GRNBoost2, seeds 1000-1099 | scripts/scenic_make_loom.py, grn_only.py; rules scenic_sketch, scenic_grn | 2026-09-08 to 09-10 | lane `looms_sketch15k/`, `adjacencies_sketch15k/` |
| 3 cisTarget | scripts/ctx_only.py; rule scenic_ctx | 2026-09-09 to 09-10 | lane `regulons_sketch15k/` |
| 4a Consensus, floor 2 | scripts/build_consensus_aucell.py; rule scenic_consensus | 2026-09-09 to 09-10 | lane `consensus_sketch15k_weighted_soft/`; its `consensus_regulons_long.tsv` was added later in September by a script that was not kept |
| 4b Consensus, floor 6 (the regulons of record) | scripts/build_consensus_aucell.py, consensus_floor6_job.sh (Slurm 52987211) | 2026-10-14 | lane `consensus_sketch15k_weighted_soft_occ6/` |
| 4c AUCell of record | scripts/cliauc_run.sh | 2026-10-26 | lane `consensus_sketch15k_weighted_soft_occ6_pyscenic_cli/` |
| 5 rho with pseudotime, Fig. 5E table | scripts/scenic_quadrant_rho_fibdec_latent.py | 2026-10-30 | `tmp/2026_10_30_5e_samap_s30fg/item1_5e/new/`; the table is staged in `xmens_final/resources/figures/panels_programs/` (sha256 pinned in `config/panels_programs.yaml`) |
| 6a Fig. 5E | scripts/plot_regulon_quadrants.py | 2026-10-30 | `xmens_figures/final_figures/Fig5_human_mouse_shared_programs/Fig5E_regulon_pseudotime_rho_human_vs_mouse_not_v9.svg` |
| 6b Fig. S28 | scripts/regulon_recurrence_cutoffs_fig.py (a Claude Science session), perrun_regulon_sizes.py | 2026-10-18 to 11-02 | `xmens_figures/final_figures/S28_scenic_supp/scenic_supp.svg` |

- `config/13_scenic.smk`: the SCENIC rules at the repository's first commit (b9fa19c, 2026-09-14), the closest kept
  version to the September run. `config/24_scenic.smk`: the same rules today (renumbered on 31 Oct 2026; since then
  they also build floor 6 and the AUCell of record).
- The 2 Nov 2026 rebuild from the final objects (`tmp/2026_11_02_final_figures/fig5e/`) reproduced the Fig. 5E table
  and the Fig. 5E SVG byte for byte; step 5 of this folder uses its version of the correlation script.
- `outputs/`: copies of the record files `compare_with_record` checks against (the Fig. 5E table, rho per regulon in
  each species, the two manuscript SVGs, and the regulon sizes of single runs behind Fig. S28).
- `make_record_tables.py` wrote `SOURCES.tsv` and `step_counts.tsv`.
