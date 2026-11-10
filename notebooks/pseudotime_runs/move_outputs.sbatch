#!/bin/bash
#SBATCH -J pt_move
#SBATCH -p eddy
#SBATCH --account=eddy_lab
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH -o /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks/pseudotime_runs/move_outputs-%j.out
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
set -euo pipefail
for f in Cevrim_XMens_Only_Fib_Dec_Cells_pseudotime.h5ad Cevrim_Human_Only_Fib_Dec_Cells_pseudotime.h5ad; do
  test ! -e /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks/$f || { echo "REFUSE: /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks/$f exists"; exit 3; }
  cp -p /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_pseudotime_notebooks/$f /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks/$f
  a=$(sha256sum /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_pseudotime_notebooks/$f | cut -d' ' -f1); b=$(sha256sum /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks/$f | cut -d' ' -f1)
  echo "$f $a $b"
  [ "$a" = "$b" ] || { echo "MISMATCH $f"; exit 4; }
done
for f in Cevrim_XMens_Only_Fib_Dec_Cells_pseudotime.h5ad Cevrim_Human_Only_Fib_Dec_Cells_pseudotime.h5ad; do rm /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_pseudotime_notebooks/$f; done
ls -la /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/notebooks /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_pseudotime_notebooks
echo MOVE_DONE
