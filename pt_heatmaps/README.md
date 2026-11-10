# xmens_pt_heatmaps: final objects -> Monocle 2 along pseudotime -> Fig. 4H, Fig. 5F and the big MI heatmaps

Three steps, no Snakemake, after a one-time install of Monocle 2: the inputs (Python), Monocle (R), and the figures,
drawn by the notebook `xmens_pt_heatmaps.ipynb`. Paths and settings live in `config.sh`. Submit from this folder;
SLURM output goes to `logs/`.

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/xmens_pt_heatmaps
export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/$(date +%Y_%m_%d)_xmens_pt_heatmaps

sbatch 0_install_monocle.sbatch  # once: monocle 2.34.0 and dplyr 1.1.4 in their own R library  -> $MONOCLE_LIB
sbatch 1_inputs.sbatch           # counts, pseudotime and the genes drawn, per species           -> $OUT/inputs/
sbatch 2_monocle.sbatch          # Monocle 2: size factors, dispersions, curves, z per gene     -> $OUT/z/
sbatch 3_figures.sbatch          # runs the notebook: Fig. 4H, Fig. 5F, the big MI heatmaps      -> $OUT/figures/
```

Each step needs the one before it finished. To queue steps 1-3 at once:

```bash
j1=$(sbatch --parsable 1_inputs.sbatch)
j2=$(sbatch --parsable --dependency=afterok:$j1 2_monocle.sbatch)
sbatch --dependency=afterok:$j2 3_figures.sbatch
```

- Step 0 is done (7 Nov 2026): `$MONOCLE_LIB` is `/n/lab_storage/eddy_lab/users/nhilgert/R/library-4.4.1-monocle2`.
  Run again, it installs only what is missing, then checks that monocle loads with dplyr 1.1.4 and runs its
  dispersion estimate on simulated data.
- Times in the checks: step 1, 1-2.5 min; step 2, 7-10 min for mouse (99 GB) and 9-14 min for human (127 GB), 21 min
  for both in one job, or as two jobs at once (`SPECIES=mouse sbatch 2_monocle.sbatch`, then `SPECIES=human ...`);
  step 3, under a minute; `compare_with_record`, half a minute.
- Steps 1-2 never overwrite: a step, or a species in step 2, whose output exists is skipped, and every output is
  written under a temporary name and renamed when complete. After a failure, submit the step again. The figures are
  drawn again each time step 3 or a notebook cell runs.
- `bash 1_inputs.sbatch --check` prints the inputs and whether they exist, without running anything.
- `xmens_pt_heatmaps.ipynb`, on the `rsc_uv` kernel, has four cells: steps 1-2 (`bash 1_inputs.sbatch && bash
  2_monocle.sbatch`, which run the R script for Monocle; skipped if `OUT` already has their outputs, otherwise they
  need a session with 180 GB), then one cell per figure, holding all of that figure's code: it draws the figure into
  `$OUT/figures/` and shows it. `OUT`, set in the first cell, is the run of record below, so the figure cells run at
  once. Step 3 runs the notebook unattended (`jupyter nbconvert --execute`) and saves the run as
  `$OUT/xmens_pt_heatmaps.executed.ipynb`.
- `figures/` keeps a copy of the run of record's outputs, because netscratch purges `$OUT` after 90 days: the SVG,
  PDF and PNG files and the MI gene-order tables, Monocle's own Fig. 4H drawing, the executed notebook and the
  comparison table, with `SHA256SUMS.txt`. They were drawn on the cluster, so their labels sit slightly differently
  from the figures of record on Drive (`xmens_figures/final_figures/`), which were drawn on a Mac (see "Differences
  from the record"). A new run writes to its own `$OUT/figures/`.
- Each script in `scripts/` also runs on its own; its header gives the usage (for step 2,
  `R_LIBS=$MONOCLE_LIB Rscript scripts/monocle_z.R SPECIES INPUTS_DIR OUT_DIR`).

## What each step does

1. `scripts/heatmap_inputs.py`: from each final object, the nuclei with a pseudotime (42,290 mouse, 47,533 human) and
   the raw counts of every gene (33,698 mouse, 36,601 human); the table of those nuclei (pseudotime, label,
   timepoint); and the genes drawn: the Yang et al. (2023) panel (`resources/yang2023_panel.tsv`, their Fig. 4H order,
   matched to mouse symbols), the 28 genes of Fig. 5F (xmens_mi step 6), and the top MI genes of each species among the
   one-to-one orthologs significant in both (1,133 mouse, 1,201 human, 373 pairs in both).
2. `scripts/monocle_z.R`, per species, everything numerical, with Monocle 2 itself:
   - `newCellDataSet()` on the raw counts, `expressionFamily = negbinomial.size()`;
   - `estimateSizeFactors()` and `estimateDispersions()` with their defaults, on all nuclei and genes (size factor =
     total counts over their geometric mean; dispersion trend a + b / mean over the genes detected in more than one
     nucleus, refitted without Cook's-distance outliers);
   - `genSmoothCurves()`, the function `plot_pseudotime_heatmap()` calls, for every gene drawn: a negative-binomial
     fit of the counts over the size factor on `sm.ns(Pseudotime, df = 3)`, with the dispersion from the trend,
     predicted at 100 evenly spaced pseudotime points. Monocle refits a gene whose fit fails with `negbinomial()`; the
     script names such genes (`fit` column: `trend` or `fallback`);
   - `plot_pseudotime_heatmap()`'s own lines, copied verbatim (its call of `genSmoothCurves()`, then log10(curve + 1)
     and z per gene; the 9 lines sit between `>>>` and `<<<` markers in `monocle_z.R`, and step 0's check finds each
     of them in the installed function, `deparse(body(monocle::plot_pseudotime_heatmap))`). Its last line, the clip at
     ±3 (its colour limits), is applied when drawing, where Monocle applies it too: the MI heatmaps order their rows
     by the unclipped z, as the record did.
   - Mouse also does Fig. 4H: the Yang panel on the 32,957 induced nuclei, pseudotime re-ranked within them (0 to 1),
     with the size factors and dispersion trend of all nuclei (Monocle keeps both when a CellDataSet is subset); and
     `plot_pseudotime_heatmap()`'s own drawing of that panel (`fig4h_monocle_plot_pseudotime_heatmap.pdf`, Monocle's
     blue-green-red colours, rows clustered), as a reference.
3. The notebook's figure cells: colours and layout only, nothing numerical. Each reads Monocle's z, clips it at ±3,
   maps it through matplotlib's `RdBu_r` and draws the panel as the figure of record has it. The two matplotlib cells
   first restore matplotlib's default settings, so a figure does not depend on the cells run before it or on
   Jupyter's display settings.
   - Fig. 4H: rows in published order; written as SVG directly: the 30 x 59.887 mm box (one rect per cell), the
     colorbar (5 pt labels) and the labelled panel (title, n, italic gene names);
   - Fig. 5F: human (left) and mouse (right), rows in drawing order;
   - big MI heatmaps: one panel per species, rows in order of pseudotime centroid (centre of mass of the positive part
     of z, ties by the point of the maximum), with the gene-order tables.

## Monocle 2

Monocle 2 is deprecated in favour of Monocle 3, but Bioconductor still ships it (monocle 2.34.0 in Bioconductor 3.20,
the release for R 4.4), and Monocle 3 has no `plot_pseudotime_heatmap()`. Its `estimateDispersions()` calls dplyr's
`group_by_()` and `select_()`, which dplyr 1.2 made defunct; step 0 therefore puts dplyr 1.1.4, in which they still run
(with deprecation warnings), in Monocle's library, which step 2 puts first on R's search path (`R_LIBS`). The lab's
library, with dplyr 1.2.1, is not changed. `$OUT/z/sessionInfo_<species>.txt` records every package version.

Comparing with Monocle's own code: `monocle_z.R` calls `newCellDataSet()`, `estimateSizeFactors()`,
`estimateDispersions()` and `genSmoothCurves()` from the installed package, and copies 9 lines of
`plot_pseudotime_heatmap()` (between its `>>>` and `<<<` markers). `sbatch 0_install_monocle.sbatch` checks those 9
lines against the installed function; its log ends `monocle_z.R: 9 lines copied from plot_pseudotime_heatmap(), each
found in the installed function`. To read the function itself, type `plot_pseudotime_heatmap` in R after
`library(monocle)` (with `R_LIBS=$MONOCLE_LIB`), or open `R/plotting.R` in the Bioconductor source,
`monocle_2.34.0.tar.gz` (release 3.20).
Cite Qiu et al. (2017), Nat Methods 14:979 (Monocle 2), and Trapnell et al. (2014), Nat Biotechnol 32:381.

## Differences from the record

`record/` holds the heatmap package of 2 Nov 2026 that drew the figures of record (its code, small inputs, z matrices,
gene-order tables and figures), with `SOURCES.tsv` (sources and sha256; each figure also matches the Drive figure byte
for byte); `record/README.md` describes how it was made. Compared with it:

- The curves and z come from Monocle 2 itself. The record used an R port of Monocle 2's code, written with VGAM,
  that does not call the monocle package, and redid the curve-to-z transform in Python for Fig. 5F and the MI
  heatmaps. The model, the settings and the transform are the same; one difference is known:
  Monocle estimates each gene's dispersion with the population variance (n), the port used the sample variance
  (n - 1). Every gene's mean matches the port's (to 10^-12), so the counts and size factors are the same; Monocle's
  empirical dispersions equal the port's converted to the population variance (to 1.3 x 10^-7), and 3 more mouse
  genes get a negative estimate, set to 0. The trend becomes a = 0.541331, b = 0.355792 for mouse (port 0.5415,
  0.3557) and a = 0.88827, b = 0.272493 for human (port 0.8883, 0.2725), and z changes by at most 1.2 x 10^-5 (below).
- Step 1 writes the raw counts of every gene for Monocle; the record wrote per-nucleus count tables of the genes drawn
  and its own dispersion table.
- The notebook's figure cells are the record's renderers cut down to what the three figures use: the drawing code is
  the same, line for line; the z transform (now Monocle's, step 2), the composition strips, the peak row order, the
  shared-gene column, the options and the provenance files are gone. The Fig. 4H colorbar with 3.6 pt labels is no
  longer drawn (the one with 5 pt labels is the figure of record).
- The notebook draws the figures on the cluster (matplotlib 3.10.8); the record drew them on a Mac (matplotlib
  3.11.1), with Helvetica Neue, which the cluster lacks. Text is stored as text in every SVG, so it reads the same,
  but every figure except the Fig. 4H box places some elements by the size of its text, measured with DejaVu Sans on
  the cluster, so positions differ slightly (the MI heatmap images come out 1 px wider), and the Fig. 5F and MI SVGs
  name the matplotlib version that drew them. Run on a Mac, the figure cells draw the figures of record (below).
- Not ported: the peak-ordered MI heatmaps (not in the paper), the PNG preview of the labelled Fig. 4H, the
  Fig. 5F z CSVs and gene table (step 2's z and step 1's gene list hold the same).

## Comparing with the record

```bash
sbatch compare_with_record.sbatch    # any time; compares the steps that have finished
```

The summary ends `logs/compare_with_record_<job ID>.out`; one row per comparison goes to
`$OUT/compare_with_record_<date>_<time>.tsv`. It compares step 1 with the package's inputs (the nuclei and their
pseudotime, labels and timepoints; the gene lists; and the counts of the Fig. 4H and Fig. 5F genes, nucleus by
nucleus), step 2's dispersion table and trend with the port's and its z matrices (clipped at ±3) with the record's
Fig. 4H and Fig. 5F z, and step 3's MI gene order and every SVG with the record's. An SVG that is not
byte-identical is compared in parts: the colours of its heatmap cells, its embedded images pixel by pixel, its text,
and the rest of the file (positions and layout).

What reproduced in the checks of November 7-8, 2026 (the run of record:
`OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_08_xmens_pt_heatmaps`, every step on the cluster, from this
folder):

- Step 1: the record's nuclei, pseudotime, labels and timepoints (42,290 mouse, 47,533 human), its gene lists (1,133
  mouse and 1,201 human MI genes, the 28 Fig. 5F genes, the Yang panel), and its counts of the Fig. 4H and Fig. 5F
  genes: 0 of 3.7 million counts differ.
- Step 2: Monocle's dispersion table has the port's 27,064 mouse and 31,735 human genes and the same means (to
  10^-12); its empirical dispersions equal the port's converted from the sample to the population variance (to
  1.3 x 10^-7). Trend: mouse a = 0.541331, b = 0.355792 (port 0.5415, 0.3557), human a = 0.88827, b = 0.272493 (port
  0.8883, 0.2725). No fallback fits, as in the record. Monocle's z differs from the record's by at most 1.2 x 10^-5
  (Fig. 4H), 6.9 x 10^-6 (Fig. 5F mouse) and 5 x 10^-6 (Fig. 5F human), and the same values are clipped.
- Step 3: the MI heatmaps have the same genes, the mouse rows in the same order, and the human rows too except RARB
  and RIPOR3 (rows 844 and 845), whose centroids differ by 10^-7 and trade places.
- Figures drawn on the cluster (step 3, running the notebook): the Fig. 4H box is byte-identical to the figure of
  record. The Fig. 4H colorbar and labelled panel and the Fig. 5F SVG have the same cell colours (the Fig. 5F
  colorbar image within 3 of 255) and the same text, at other positions (DejaVu Sans metrics, matplotlib 3.10.8). The
  MI SVG has the same text; its heatmap images are 1 px wider than the record's, so they are not compared pixel by
  pixel.
- The notebook's figure cells, run on a Mac under Jupyter's plotting backend from the same z (matplotlib 3.11.1,
  Helvetica Neue): the three Fig. 4H SVGs and the Fig. 5F SVG are byte-identical to the figures of record; the MI SVG
  differs only in the two human rows that trade places and in 392 pixels of the mouse panel (by at most 3 of 255), and
  in the image ids matplotlib derives from the pixels.
- The figure cells hold the code of the former drawing scripts (`scripts/fig4h.py`, `fig5f.py`, `mi_heatmaps.py`):
  all 11 figure files they draw are byte-identical to the scripts' on the cluster (Jupyter kernel) and on the Mac.
- `z/fig4h_monocle_plot_pseudotime_heatmap.pdf` is Monocle's own drawing of the Fig. 4H panel (its blue-green-red
  colours, rows clustered), for reference.
- Run to run (the last two runs), every input, z matrix, dispersion table and figure is byte-identical (gzip files
  compared decompressed), except Monocle's PDF, which records when it was drawn.
