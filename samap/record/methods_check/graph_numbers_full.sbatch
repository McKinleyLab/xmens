#!/bin/bash
#SBATCH -p eddy -A eddy_lab -c 2 --mem=48G -t 0:30:00 -J samap_methods_check_full
#SBATCH -o /n/lab_storage/eddy_lab/users/nhilgert/tmp/2026_11_06_xmens_samap_build/methods_checks/%x_%j.out
#SBATCH --exclude=holy8a28510,holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
set -euo pipefail
export HDF5_USE_FILE_LOCKING=FALSE
/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python /n/lab_storage/eddy_lab/users/nhilgert/tmp/2026_11_06_xmens_samap_build/methods_checks/graph_numbers_full.py /n/lab_storage/eddy_lab/users/nhilgert/tmp/2026_11_06_xmens_samap_build/methods_checks/graph_numbers_full.json

