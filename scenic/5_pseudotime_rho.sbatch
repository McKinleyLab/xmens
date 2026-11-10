#!/bin/bash
# Step 5 of 6: regulon activity vs pseudotime, both species, in one task.
#   scripts/regulon_rho_pseudotime.py reads the activity of record (step 4: consensus_..._occ6_pyscenic_cli/
#   auc_mtx_weighted.csv) and obs['pseudotime'] of the final objects ($FINAL_MOUSE, $FINAL_HUMAN in config.sh), leaves
#   out nuclei without pseudotime (the 907 apoptotic mouse decidual nuclei), and writes $OUT/correlations/:
#     regulon_rho_pseudotime_mouse.tsv, regulon_rho_pseudotime_human.tsv
#                                     Spearman rho of every regulon's activity with pseudotime (193 and 141 regulons)
#     regulon_quadrants_input.csv     the 111 TFs with a regulon in both species (Ensembl 116 one-to-one pairs):
#                                     human_tf, mouse_tf, rho_human, rho_mouse, label (1 = one of the 20 canonical
#                                     decidualization TFs); the table Fig. 5E is drawn from
#     provenance.json                 sha256 of every input, and checks (nuclei per species, pseudotime per subtype)
#   Logs: $OUT/logs/5_pseudotime_rho.log, logs/5_pseudotime_rho_<job ID>.out
#   Settings: every setting of each Python call, library defaults included, in $OUT/logs/<log name>.settings.json
#     (scripts/run_with_settings.py; see README.md)
# The folder is written as correlations.partial.<job ID> and renamed when complete; if it exists, nothing is done.
#
# Usage:
#   bash 5_pseudotime_rho.sbatch --check     inputs and output found; run nothing
#   sbatch 5_pseudotime_rho.sbatch

#SBATCH --job-name=5_pseudotime_rho
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
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

# config.sh sets $OUT, $PY, $FINAL_MOUSE, $FINAL_HUMAN, $ORTHOLOG_PAIRS, $AUCELL and $LABELLED.
source "$HERE/config.sh"

ACT_MOUSE=$OUT/scenic_seeds/mouse/$AUCELL/auc_mtx_weighted.csv
ACT_HUMAN=$OUT/scenic_seeds/human/$AUCELL/auc_mtx_weighted.csv
DEST=$OUT/correlations

if [ "${1:-}" = "--check" ]; then
    echo "activity: mouse $(exists "$ACT_MOUSE"), human $(exists "$ACT_HUMAN")"
    echo "final objects: mouse $(exists "$FINAL_MOUSE"), human $(exists "$FINAL_HUMAN"); orthologs $(exists "$ORTHOLOG_PAIRS")"
    echo "output $DEST: $(exists "$DEST")"
    exit 0
fi

if [ -e "$DEST" ]; then
    echo "$DEST exists; nothing to do"
    exit 0
fi
TMP=$DEST.partial.${SLURM_JOB_ID:-$$}
mkdir -p "$OUT/logs"
echo "started $(date) on $(hostname)"

# --expected-n is recorded in provenance.json, not enforced: the nuclei with a pseudotime value in the final objects.
"$PY" "$HERE/scripts/run_with_settings.py" \
    --settings "$OUT/logs/5_pseudotime_rho.settings.json" \
    -- "$HERE/scripts/regulon_rho_pseudotime.py" \
    --activity-mouse "$ACT_MOUSE" \
    --activity-human "$ACT_HUMAN" \
    --object-mouse "$FINAL_MOUSE" \
    --object-human "$FINAL_HUMAN" \
    --pt-col pseudotime \
    --orthologs "$ORTHOLOG_PAIRS" \
    --labelled "$LABELLED" \
    --expected-n mouse=42290,human=47533 \
    --out-dir "$TMP" \
    > "$OUT/logs/5_pseudotime_rho.log" 2>&1
mv "$TMP" "$DEST"

# The script's summary line (TF pairs, labelled TFs) also goes to the SLURM log.
tail -n 1 "$OUT/logs/5_pseudotime_rho.log"
echo "finished $(date)"
