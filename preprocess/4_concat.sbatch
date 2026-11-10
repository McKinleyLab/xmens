#!/bin/bash
# Step 4 of 5: per-sample QC'd objects -> one object per species.
#   Reads:  $OUT/qc/<sample>.h5ad (step 3)
#   Writes: $OUT/qc_passed_mouse.h5ad and $OUT/qc_passed_human.h5ad
#   Logs:   $OUT/logs/4_concat_<species>.log and logs/4_concat_<job ID>.out
# A species whose object exists is skipped, and a missing per-sample object stops the job.
#
# Usage:
#   sbatch 4_concat.sbatch

#SBATCH --job-name=4_concat
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out  # job name and job ID; relative to the folder sbatch was run from

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# The pipeline folder holds config.sh and samples.tsv. Run with bash, this script sits in it; under sbatch, SLURM
# runs a copy from its spool folder, and the pipeline folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the preprocess folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh sets $OUT, $PY, N_SAMPLES, and field <row> <column>, which prints one value from samples.tsv.
source "$HERE/config.sh"

for species in mouse human; do
    out=$OUT/qc_passed_$species.h5ad
    if [ -e "$out" ]; then
        echo "$species: $out exists; skipped"
        continue
    fi

    # The species' per-sample objects in samples.tsv order, which sets the order of nuclei in the object, as in
    # the run of record.
    inputs=()
    for row in $(seq 1 "$N_SAMPLES"); do
        if [ "$(field "$row" species)" = "$species" ]; then
            inputs+=("$OUT/qc/$(field "$row" sample_id).h5ad")
        fi
    done

    "$PY" "$HERE/scripts/concat.py" \
        --species "$species" \
        --out "$out" \
        --log "$OUT/logs/4_concat_$species.log" \
        "${inputs[@]}"

    # The log's last line, the file written, also goes to the SLURM log.
    tail -n 1 "$OUT/logs/4_concat_$species.log"
done
