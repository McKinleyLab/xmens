#!/bin/bash
# Step 6 of 6: Fig. 5E and Fig. S28, in one task.
#   Fig. 5E: scripts/plot_regulon_quadrants.py draws step 5's regulon_quadrants_input.csv: rho in human (x) and mouse
#     (y) for the 111 TF pairs, coloured when |rho| >= 0.2 in both species, the 13 TFs of the manuscript panel named.
#     -> $OUT/figures/Fig5E_regulon_quadrants.{svg,pdf,png,provenance.json}
#   Fig. S28: scripts/s28_consensus_edges.py lists the links of step 4's floor-2 consensus (one row per TF-target pair,
#     with the number of runs containing it) and scripts/s28_perrun_sizes.py the size of every regulon in each of the
#     100 runs (from $REGULONS_FROM), into $OUT/figures/s28_inputs/; scripts/plot_regulon_recurrence_cutoffs.py then
#     draws the TF and target thresholds -> $OUT/figures/FigS28_regulon_recurrence_cutoffs.{svg,png,provenance.json}
#   The two regulon tables run in the pySCENIC environment, because the regulons are pickled pySCENIC objects; the
#   plots run in rsc_uv. The scripts ask for Helvetica Neue, then Helvetica, Arial, Nimbus Sans; on the cluster they
#   get Nimbus Sans.
#   Logs: $OUT/logs/6_fig5e.log, 6_figS28.log, logs/6_figures_<job ID>.out
#   Settings: every setting of each Python call, library defaults included, in $OUT/logs/<log name>.settings.json
#     (scripts/run_with_settings.py; see README.md)
# The folder is written as figures.partial.<job ID> and renamed when complete; if it exists, nothing is done.
#
# Usage:
#   bash 6_figures.sbatch --check     inputs and output found; run nothing
#   sbatch 6_figures.sbatch

#SBATCH --job-name=6_figures
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G  # the regulons of all 200 runs are read one at a time
#SBATCH --time=01:00:00
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506  # GPU nodes
#SBATCH --output=logs/%x_%j.out  # job name and job ID; relative to the folder sbatch was run from

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# The pipeline folder holds config.sh and runs.tsv. Run with bash, this script sits in it; under sbatch, SLURM runs a
# copy from its spool folder, and the pipeline folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the scenic folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, $REGULONS_FROM, $PY, $PY_SCENIC, $ORTHOLOG_COMPONENTS, $TAG, $CONS6, $CONS2, $NAMED and
# $RHO_THRESHOLD.
source "$HERE/config.sh"

TABLE=$OUT/correlations/regulon_quadrants_input.csv
DEST=$OUT/figures

if [ "${1:-}" = "--check" ]; then
    echo "Fig. 5E table $(exists "$TABLE")"
    for sp in mouse human; do
        n_pkl=$(ls "$REGULONS_FROM/scenic_seeds/$sp/regulons_$TAG"/regulons_seed*.pkl 2>/dev/null | wc -l | tr -d ' ' || true)
        echo "$sp: floor-6 report $(exists "$OUT/scenic_seeds/$sp/$CONS6/recurrence_report.tsv")," \
            "floor-2 regulons $(exists "$OUT/scenic_seeds/$sp/$CONS2/consensus_regulons.pkl"), $n_pkl run regulon files"
    done
    echo "output $DEST: $(exists "$DEST")"
    exit 0
fi

if [ -e "$DEST" ]; then
    echo "$DEST exists; nothing to do"
    exit 0
fi
TMP=$DEST.partial.${SLURM_JOB_ID:-$$}
mkdir -p "$TMP/s28_inputs" "$OUT/logs"
echo "started $(date) on $(hostname)"

# Fig. 5E. The plot script checks the table against a sha256 and records where it came from; here that is step 5's
# table itself (compare_with_record.sbatch checks it against the table of record).
"$PY" "$HERE/scripts/run_with_settings.py" \
    --settings "$OUT/logs/6_fig5e.settings.json" \
    -- "$HERE/scripts/plot_regulon_quadrants.py" \
    --table "$TABLE" \
    --table-sha256 "$(sha256sum "$TABLE" | cut -c1-64)" \
    --table-provenance "$OUT/correlations/provenance.json" \
    --source-of-record "$TABLE" \
    --threshold "$RHO_THRESHOLD" \
    --named "$NAMED" \
    --out "$TMP/Fig5E_regulon_quadrants" \
    > "$OUT/logs/6_fig5e.log" 2>&1

# Fig. S28: the floor-2 consensus links of both species, the regulon sizes of the 100 runs, then the figure.
LOG=$OUT/logs/6_figS28.log
for sp in mouse human; do
    "$PY_SCENIC" "$HERE/scripts/run_with_settings.py" \
        --settings "$OUT/logs/6_figS28_links_$sp.settings.json" \
        -- "$HERE/scripts/s28_consensus_edges.py" \
        --pkl "$OUT/scenic_seeds/$sp/$CONS2/consensus_regulons.pkl" \
        --out "$TMP/s28_inputs/${sp}_consensus_regulons_long.tsv" \
        >> "$LOG" 2>&1
done
"$PY_SCENIC" "$HERE/scripts/run_with_settings.py" \
    --settings "$OUT/logs/6_figS28_run_sizes.settings.json" \
    -- "$HERE/scripts/s28_perrun_sizes.py" \
    --lane "$REGULONS_FROM/scenic_seeds" \
    --regulon-dir "regulons_$TAG" \
    --n-expected 100 \
    --out-json "$TMP/s28_inputs/perrun_sizes.json" \
    >> "$LOG" 2>&1
"$PY" "$HERE/scripts/run_with_settings.py" \
    --settings "$OUT/logs/6_figS28.settings.json" \
    -- "$HERE/scripts/plot_regulon_recurrence_cutoffs.py" \
    --report-mouse "$OUT/scenic_seeds/mouse/$CONS6/recurrence_report.tsv" \
    --report-human "$OUT/scenic_seeds/human/$CONS6/recurrence_report.tsv" \
    --long-mouse "$TMP/s28_inputs/mouse_consensus_regulons_long.tsv" \
    --long-human "$TMP/s28_inputs/human_consensus_regulons_long.tsv" \
    --perrun "$TMP/s28_inputs/perrun_sizes.json" \
    --orthologs "$ORTHOLOG_COMPONENTS" \
    --out "$TMP/FigS28_regulon_recurrence_cutoffs" \
    >> "$LOG" 2>&1
mv "$TMP" "$DEST"
ls -l "$DEST"
echo "finished $(date)"
