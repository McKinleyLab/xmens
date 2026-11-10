#!/bin/bash
# Compare md5 checksums of the Cevrim h5ad backup originals against the netscratch copies.
# Read-only for all h5ad files; the only file written is the log below.

DATA=/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data
SCRATCH=/n/netscratch/eddy_lab/Everyone/nhilgert
LOG=/n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/checksums/cevrim_checksums_$(date +%Y%m%d).txt

{
  echo "# Cevrim h5ad checksum comparison"
  echo "# date: $(date)   host: $(hostname)"
  echo "# originals: $DATA"
  echo "# copies:    $SCRATCH"
  echo
  ls -l "$DATA"/20251008_Cevrim*.h5ad "$SCRATCH"/20251008_Cevrim*.h5ad | sed 's/^/# /'
  echo
  for path in "$DATA"/20251008_Cevrim*.h5ad; do
    f=$(basename "$path")
    a=$(md5sum "$DATA/$f")
    b=$(md5sum "$SCRATCH/$f")
    echo "$a"
    echo "$b"
    if [ -z "$a" ] || [ -z "$b" ]; then
      echo "# ERROR (could not read one or both): $f"
    elif [ "${a%% *}" = "${b%% *}" ]; then
      echo "# MATCH: $f"
    else
      echo "# DIFFER: $f"
    fi
    echo
  done
} 2>&1 | tee "$LOG"
