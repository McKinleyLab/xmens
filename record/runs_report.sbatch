#!/bin/bash
# Fills in the report blocks of <component>/record/RUNS.md for the runs in <component>/record/jobs.tsv, once they have
# ended: each step's sacct record, a sha256 manifest of every file in each run folder, and each rerun compared with
# its run file by file (record/runs_report.py has the details).
#   Writes: <component>/record/RUNS.md (the report blocks only), <component>/record/manifests/ and reports/
#   Log:    record/runs_report_<job ID>.out
#
# Usage, from the repository root, queued after the last job of every chain it reports on:
#   sbatch --dependency=afterany:<job>:<job>... record/runs_report.sbatch <component>... [--reference RUN=FOLDER]
#SBATCH --job-name=runs_report
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G  # compares h5ad files one dataset at a time; the largest datasets are a few GB
#SBATCH --time=04:00:00  # hashing about 50 GB of run folders, and comparing the files whose bytes differ
#SBATCH --output=record/runs_report_%j.out  # relative to the folder sbatch was run from (the repository root)

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
export PYTHONDONTWRITEBYTECODE=1  # no __pycache__ in the repository
PY=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python

echo "runs_report on $(hostname), $(date -Iseconds), repository at $(git rev-parse --short HEAD)"
"$PY" record/runs_report.py "$@"
echo "done $(date -Iseconds)"
