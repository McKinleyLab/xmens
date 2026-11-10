# The 2026 Snakemake pipeline

The analysis for the revision was developed from June to November 2026 as a Snakemake pipeline,
`xmens_final`. It lives in its own repository with its complete commit history:

https://github.com/nhilgert/xmens_final (tag `snakemake-final`, 2026-11-09)

Starting in November 2026 the analysis moves into this repository one step at a time, as plain
sbatch scripts ("one script per step"). Each step folder carries a `record/SOURCES.tsv` naming the
`xmens_final` commit and script path it was copied from, plus the run of record and checksums, so
every result can be traced back to the pipeline commit that first produced it.

The October 2026 repository `xmens-pipeline` was a simplified snapshot of `xmens_final` and is
archived; `xmens_final` is the reference for the pipeline's history.
