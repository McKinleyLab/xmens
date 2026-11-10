# Where the run of record's logs and outputs are

Lane = `/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/scenic_v2_2026_08_28`; `NH` =
`/n/lab_storage/eddy_lab/users/nhilgert`.

## Steps 1-3 and the floor-2 consensus (Snakemake, 7-10 Sep 2026)

- Rule logs: `NH/xmens_final/logs/scenic/`: `pooled_<species>.log`, `sketch_<species>_seed<seed>.log`,
  `grn_<species>_seed<seed>.log`, `ctx_<species>_seed<seed>.log` (200 each), `consensus_<species>.log`.
- Snakemake's own SLURM logs: lane `slurm_logs/` (`smk_scenic_all_<job ID>.{out,err}`).
- Outputs: lane `scenic_seeds/<species>/`: `<species>_dsc_pooled.{h5ad,loom,manifest.json}`, `looms_sketch15k/`
  (a loom and a manifest per run), `adjacencies_sketch15k/`, `regulons_sketch15k/`,
  `consensus_sketch15k_weighted_soft/`.
- The lane's `README.md` records the decisions of September and October (15k sketches, the consensus variants,
  floor 6 as the primary consensus on 15 Oct).

## Floor 6, 14 Oct 2026

- Job script: `scripts/consensus_floor6_job.sh` here (Slurm 52987211; 16.5 min, 48 GB for both species).
- Logs: lane `consensus_variant_logs/consensus_occ6_{mouse,human}_20261014.log`.
- Outputs: lane `scenic_seeds/<species>/consensus_sketch15k_weighted_soft_occ6/`.

## AUCell of record, 26 Oct 2026

- Job script: `scripts/cliauc_run.sh` here.
- Outputs and a README with the exact commands: lane
  `scenic_seeds/<species>/consensus_sketch15k_weighted_soft_occ6_pyscenic_cli/` (`consensus_regulons.yaml`,
  `auc_mtx_weighted.csv` = the activity of record, `auc_mtx_unweighted.csv`).

## Fig. 5E table and panel, 30 Oct 2026; rebuilt 2 Nov 2026

- 30 Oct: `NH/tmp/2026_10_30_5e_samap_s30fg/item1_5e/new/` (`rho_fibdec_latent.log`, the provenance JSON, the
  script and its diff against the repository version); the panel script `plot_regulon_quadrants.py` one folder up.
- 2 Nov rebuild from the final objects: `NH/tmp/2026_11_02_final_figures/fig5e/` (`README.md`, `checks/`,
  `out/regulon_table/`, `out/panel_5E_one.svg`), byte-identical to the 30 Oct table and panel.
- Manuscript panel: `NH/xmens_figures/final_figures/Fig5_human_mouse_shared_programs/`.

## Fig. S28, 18 Oct - 2 Nov 2026

- Drawn in a Claude Science session from the floor-6 recurrence reports, the floor-2 link tables and the per-run
  regulon sizes (`scripts/regulon_recurrence_cutoffs_fig.py`, `scripts/perrun_regulon_sizes.py` here).
- Manuscript file: `NH/xmens_figures/final_figures/S28_scenic_supp/scenic_supp.svg`.
- The repository's port of the script and its outputs: `NH/xmens_final/results/figures/panels/scenic_params/`
  (same text and numbers as the manuscript file).
