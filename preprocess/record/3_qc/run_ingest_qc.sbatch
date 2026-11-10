#!/bin/bash
#SBATCH -p eddy
#SBATCH --account=eddy_lab
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o slurm-qc-rerun-%j.out
# Ingest -> QC -> concat with the pipeline's own scripts, outside Snakemake (run_ingest_qc.py), then the checks.
# Submitted by verify_qc.sh; three modes:
#   SAMPLE set           ingest + QC of that one sample (verify_qc.sh submits one job per sample)
#   array task           ingest + QC of the sample on data row (task + 1) of samples.tsv
#   neither              all samples (finished ones are skipped), concat per species, then three checks in $OUT:
#                          compare_<stamp>.txt            vs the pipeline's objects (compare_qc_outputs.py)
#                          vs_2025_qc_<stamp>/            vs the 2025 per-lane QC (compare_qc_vs_2025.py)
#                          vs_2025_final_<species>_<stamp>.txt   vs the 2025 final objects (check_qc_vs_submission.py)
# INPUTS=andy: the pipeline's own inputs (xmens_inputs); INPUTS=regen: the February 2026 CellBender re-run.
set -euo pipefail
INPUTS=${INPUTS:?set INPUTS=andy or INPUTS=regen}
PY=${PY:-/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python}
OUT=${OUT:-/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_05_qc_rerun_${INPUTS}}
SAMPLES=/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/config/samples.tsv
export HDF5_USE_FILE_LOCKING=FALSE
cd "${SLURM_SUBMIT_DIR:-.}"
echo "node $(hostname); CPUs allocated ${SLURM_CPUS_PER_TASK:-?} of $(getconf _NPROCESSORS_ONLN)"
if [[ -z "${SAMPLE:-}" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    SAMPLE=$(awk -F'\t' -v row="$((SLURM_ARRAY_TASK_ID + 2))" 'NR == row {print $1}' "$SAMPLES")
    [[ -n "$SAMPLE" ]] || { echo "no sample for array task $SLURM_ARRAY_TASK_ID"; exit 2; }
fi
if [[ -n "${SAMPLE:-}" ]]; then
    "$PY" run_ingest_qc.py --inputs "$INPUTS" --out "$OUT" --samples "$SAMPLE" --steps ingest qc
else
    "$PY" run_ingest_qc.py --inputs "$INPUTS" --out "$OUT" --steps ingest qc concat
    STAMP=$(date +%Y%m%d_%H%M%S)
    status=0   # run every check even if one reports a FAIL; the job's exit status says whether any did
    "$PY" compare_qc_outputs.py --new "$OUT" | tee "$OUT/compare_$STAMP.txt" || status=1
    "$PY" compare_qc_vs_2025.py --inputs "$INPUTS" --current "$OUT/qc" --out "$OUT/vs_2025_qc_$STAMP" || status=1
    for species in mouse human; do
        "$PY" check_qc_vs_submission.py --species "$species" --inputs "$INPUTS" --current "$OUT/qc" \
            --out "$OUT/vs_2025_final_${species}_lost_$STAMP.tsv" | tee "$OUT/vs_2025_final_${species}_$STAMP.txt" \
            || status=1
    done
    exit $status
fi
