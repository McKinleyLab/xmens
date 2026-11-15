# pt_heatmaps: expression along pseudotime with Monocle 2 (Fig. 4H, Fig. 5F, the MI heatmaps)

Smoothed expression curves along pseudotime from Monocle 2's `genSmoothCurves()`, z-scored per gene as its
`plot_pseudotime_heatmap()` does, and the heatmaps drawn from them by `xmens_pt_heatmaps.ipynb`. Paths and settings
are in `config.sh`.

## Run

```bash
cd pt_heatmaps && mkdir -p logs
export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_pt_heatmaps    # a new folder for each run

sbatch 0_install_monocle.sbatch                                      # once: monocle 2.34.0 in its own R library
j1=$(sbatch --parsable 1_inputs.sbatch)                              # counts, pseudotime and the genes drawn
j2=$(sbatch --parsable --dependency=afterok:$j1 2_monocle.sbatch)    # Monocle 2: curves and z per gene
sbatch --dependency=afterok:$j2 3_figures.sbatch                     # runs the notebook, which draws the figures
```

Outputs go to `$OUT/inputs/`, `$OUT/z/` and `$OUT/figures/`. A step whose output exists is skipped.

## Notes

- Step 2 follows Monocle 2's defaults: `estimateSizeFactors()` and `estimateDispersions()` on all nuclei and genes,
  then `genSmoothCurves()` (negative binomial, `sm.ns(Pseudotime, df = 3)`, 100 points), log10(curve + 1) and z per
  gene. z is clipped at +/-3 when drawn.
- Fig. 4H uses the induced mouse nuclei only, with pseudotime re-ranked within them.
- Monocle 2's `estimateDispersions()` needs dplyr 1.1.4, so step 0 installs both into a separate R library
  (`$MONOCLE_LIB`). `$OUT/z/sessionInfo_<species>.txt` records every package version.
- `figures/` keeps the executed notebook and the MI heatmaps' gene order from an earlier run (8 Nov 2026).

## Record

`record/` holds the earlier heatmap code and its outputs, described in `record/README.md`.
