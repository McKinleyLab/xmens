#!/bin/bash
#SBATCH --job-name=phmmer_to_symbols
#SBATCH --partition=eddy
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/n/eddy_lab/users/nhilgert/decidualization/revisions/samap_integrate/logs/phmmer_to_symbols_%j.out
#SBATCH --error=/n/eddy_lab/users/nhilgert/decidualization/revisions/samap_integrate/logs/phmmer_to_symbols_%j.err

set -euo pipefail

WORKDIR=/n/eddy_lab/users/nhilgert/decidualization/revisions/samap_integrate
cd "${WORKDIR}"
mkdir -p logs homology/hsmm_hmmer

# --- Environment ---
source /n/eddy_lab/users/nhilgert/.uv/envs/saturn_uv/bin/activate

echo "=== JOB INFO ==="
echo "Job ID:    ${SLURM_JOB_ID}"
echo "Host:      $(hostname)"
echo "Working:   $(pwd)"
echo "Python:    $(which python)"
echo "Start:     $(date)"
echo ""

echo "=== Direction 1: hu_vs_mu (human queries vs mouse targets) ==="
time python scripts/phmmer_to_samap_symbols.py hu_vs_mu
echo ""

echo "=== Direction 2: mu_vs_hu (mouse queries vs human targets) ==="
time python scripts/phmmer_to_samap_symbols.py mu_vs_hu
echo ""

echo "=== OUTPUT FILES ==="
ls -lh homology/hsmm_hmmer/

echo ""
echo "=== SANITY: known orthologs in hs_to_mm.txt ==="
for pair in "HAND2 Hand2" "FOXO1 Foxo1" "IL11 Il11" "WNT5A Wnt5a" "PGR Pgr"; do
    set -- $pair
    line=$(grep -P "^$1\t$2\b" homology/hsmm_hmmer/hs_to_mm.txt || true)
    if [[ -n "$line" ]]; then
        echo "  ✓ $1 ↔ $2  →  $line"
    else
        echo "  ✗ $1 ↔ $2  NOT FOUND"
    fi
done

echo ""
echo "=== DONE ==="
echo "End:       $(date)"
