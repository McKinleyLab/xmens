#!/bin/bash
# MD5 checksums of the paper's raw FASTQs: the 15 gene-expression libraries (samples.tsv, column fastq_dir) and the
# 15 spatial-barcode libraries (spatial_libraries.tsv, column spatial_fastq_dir), one line per file, as GEO asks for
# them. Reads every file once; changes nothing under $FASTQS.
#   Writes: record/fastq_md5_<date>.tsv with the columns
#           sample_id  lane  library (gex|spatial)  fastq_dir  file  bytes  md5
#           A file of that name that already exists is never overwritten (the job stops instead).
#   Log:    logs/checksum_fastqs_<job ID>.out
# About 1 TB of FASTQs, 8 files at a time: expect a few hours.
#
# Usage:
#   bash checksum_fastqs.sbatch --check     list the libraries, folders and file counts; hash nothing
#   sbatch checksum_fastqs.sbatch           hash everything
#   bash checksum_fastqs.sbatch             the same, in the current shell (an interactive session, not a login node)

#SBATCH --job-name=checksum_fastqs
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=8  # md5sum runs 8 files at a time (THREADS below)
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out  # job name and job ID; relative to the folder sbatch was run from

# Stop at the first failed command, unset variable or failed stage of a pipe.
set -euo pipefail

# The pipeline folder holds config.sh and the two sample sheets. Run with bash, this script sits in it; under sbatch,
# SLURM runs a copy from its spool folder, and the pipeline folder is then the one sbatch was run from.
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -e "$HERE/config.sh" ]; then
    HERE=${SLURM_SUBMIT_DIR:-$PWD}
fi
if [ ! -e "$HERE/config.sh" ]; then
    echo "submit this from the preprocess folder (no config.sh in $HERE)" >&2
    exit 2
fi

# config.sh gives $FASTQS and field <row> <column>. It insists on $OUT, the run folder, which this script does not
# use (it writes only under $HERE/record), so a placeholder satisfies it.
export OUT=${OUT:-/dev/null}
source "$HERE/config.sh"

THREADS=${SLURM_CPUS_PER_TASK:-8}
SPATIAL=$HERE/spatial_libraries.tsv
if [ ! -e "$SPATIAL" ]; then
    echo "no $SPATIAL (sample_id, lane, spatial_fastq_dir, puck_id)" >&2
    exit 2
fi

# spatial_dir <lane>: the lane's spatial FASTQ folder from spatial_libraries.tsv (same lookup as field, other file).
spatial_dir() {
    awk -F'\t' -v lane="$1" 'NR > 1 && $2 == lane { print $3 }' "$SPATIAL"
}

# One line per library: sample_id, lane, library, fastq_dir. Both sheets are keyed by lane; a lane missing from
# spatial_libraries.tsv stops the script.
libraries() {
    for row in $(seq 1 "$N_SAMPLES"); do
        sample=$(field "$row" sample_id)
        lane=$(field "$row" lane)
        printf '%s\t%s\tgex\t%s\n' "$sample" "$lane" "$(field "$row" fastq_dir)"
        spatial=$(spatial_dir "$lane")
        if [ -z "$spatial" ]; then
            echo "$lane: not in $SPATIAL" >&2
            exit 1
        fi
        printf '%s\t%s\tspatial\t%s\n' "$sample" "$lane" "$spatial"
    done
}

# --check: the libraries, their folders and how many FASTQs each holds; then stop.
if [ "${1:-}" = "--check" ]; then
    libraries | while IFS=$'\t' read -r sample lane library dir; do
        n=$(ls "$FASTQS/$dir"/*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
        echo "$sample $library: $FASTQS/$dir, $n FASTQ files"
    done
    exit 0
fi

OUT_TSV=$HERE/record/fastq_md5_$(date +%Y%m%d).tsv
if [ -e "$OUT_TSV" ]; then
    echo "$OUT_TSV exists; not overwriting it" >&2
    exit 1
fi
mkdir -p "$HERE/record"

# Work files next to the output (never in the home directory); removed on exit.
WORK=$(mktemp -d "$HERE/record/.checksum_fastqs.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# files.tsv: sample_id, lane, library, fastq_dir, full path, for every FASTQ of every library.
libraries | while IFS=$'\t' read -r sample lane library dir; do
    for path in "$FASTQS/$dir"/*.fastq.gz; do
        printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "$lane" "$library" "$dir" "$path"
    done
done > "$WORK/files.tsv"
echo "$(wc -l < "$WORK/files.tsv") files in $(( N_SAMPLES * 2 )) libraries; hashing $THREADS at a time, started $(date)"

# md5sum prints "<md5>  <path>"; the paths come back in completion order, so they are joined to files.tsv by path.
cut -f5 "$WORK/files.tsv" | xargs -P "$THREADS" -n 1 md5sum > "$WORK/md5.txt"

{
    printf 'sample_id\tlane\tlibrary\tfastq_dir\tfile\tbytes\tmd5\n'
    awk -F'\t' 'NR == FNR { split($0, parts, "  "); md5[parts[2]] = parts[1]; next }
                { path = $5; n = split(path, bits, "/"); file = bits[n]
                  "stat -c %s " path | getline bytes; close("stat -c %s " path)
                  if (!(path in md5)) { print "no md5 for " path > "/dev/stderr"; exit 1 }
                  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, file, bytes, md5[path] }' \
        "$WORK/md5.txt" "$WORK/files.tsv"
} > "$OUT_TSV"

echo "wrote $OUT_TSV: $(( $(wc -l < "$OUT_TSV") - 1 )) files, finished $(date)"
