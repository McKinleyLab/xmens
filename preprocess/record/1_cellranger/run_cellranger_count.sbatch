#!/bin/bash
#SBATCH --job-name=cellranger_count
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-15
#SBATCH --output=/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_cellranger_count/logs/%x_%A_%a.out
#SBATCH --error=/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_cellranger_count/logs/%x_%A_%a.err

# Cell Ranger 7.2.0 count matrices for the 15 X-Mens lanes, from the FASTQs in xmens_final/dataverse/fastqs.
# References and settings are those of the February 2026 re-run (revisions/velocity/<lane>/cellranger/_cmdline):
# transgenic_mus (mm10 + GsD + GqD) for mouse, GRCh38-2024-A for human, 8 cores, 96 GB, chemistry auto,
# introns included (the 7.2.0 default). --no-bam skips the BAM, which does not change the matrices.
#
#   bash run_cellranger_count.sbatch --check    # print each lane's command and check its FASTQs; runs nothing
#   sbatch run_cellranger_count.sbatch          # one array task per lane
#
# Each lane writes only into OUT/<lane>. A lane whose outs/raw_feature_bc_matrix.h5 already exists is skipped;
# rerunning a failed lane resumes it in place.
set -euo pipefail

CELLRANGER=/n/lab_storage/eddy_lab/users/nhilgert/software/cellranger-7.2.0/cellranger
FASTQS=/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/dataverse/fastqs
MOUSE_REF=/n/lab_storage/eddy_lab/users/nhilgert/reference/transgenic_mus
HUMAN_REF=/n/lab_storage/eddy_lab/users/nhilgert/reference/human/refdata-gex-GRCh38-2024-A
OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_cellranger_count
CACHE=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch/cellranger

# lane, FASTQ folder under FASTQS, sample index (FASTQ name prefix), species, reads in the February run
LANES="mck_8a  FG/SI-TT-F7    SI-TT-F7   mouse  837,511,274
mck_8b  FG/SI-TT-F8    SI-TT-F8   mouse  827,006,620
mck_1   mck_1/SI-TT-A5 SI-TT-A5   mouse  191,981,288
mck_12a FG/SI-TT-F4    SI-TT-F4   mouse  586,181,794
mck_12b FG/SI-TT-F5    SI-TT-F5   mouse  364,217,847
mck_15  FG/SI-TT-G6    SI-TT-G6   mouse  230,403,475
mck_10a FG/SI-TT-G4    SI-TT-G4   mouse  384,897,848
mck_10b FG/SI-TT-G5    SI-TT-G5   mouse  790,967,709
mck_11a FG/SI-TT-F1    SI-TT-F1   mouse  1,039,108,587
mck_6a  AE/SI-TT-A8    SI-TT-A8   human  992,622,513
mck_6b  AE/SI-TT-A9    SI-TT-A9   human  894,677,488
mck_6c  AE/SI-TT-A10   SI-TT-A10  human  790,009,766
mck_5a  AE/SI-TT-A5    SI-TT-A5   human  763,434,770
mck_5b  AE/SI-TT-A6    SI-TT-A6   human  766,666,696
mck_5c  AE/SI-TT-A7    SI-TT-A7   human  703,360,874"

count_command() {   # the cellranger count call for line $1 of LANES
    read -r LANE DIR SAMPLE SPECIES READS <<< "$(sed -n "${1}p" <<< "$LANES")"
    REF=$MOUSE_REF
    if [ "$SPECIES" = human ]; then REF=$HUMAN_REF; fi
    COMMAND=("$CELLRANGER" count --id="$LANE" --transcriptome="$REF" --fastqs="$FASTQS/$DIR" --sample="$SAMPLE"
             --output-dir="$OUT/$LANE" --localcores=8 --localmem=96 --no-bam --disable-ui)
}

if [ "${1:-}" = "--check" ]; then
    for i in $(seq 1 15); do
        count_command "$i"
        r1=$(ls "$FASTQS/$DIR/${SAMPLE}"_S*_L00*_R1_001.fastq.gz 2>/dev/null | wc -l || true)
        r2=$(ls "$FASTQS/$DIR/${SAMPLE}"_S*_L00*_R2_001.fastq.gz 2>/dev/null | wc -l || true)
        state="to run"
        if [ -e "$OUT/$LANE/outs/raw_feature_bc_matrix.h5" ]; then state="done, skipped"; fi
        echo "$LANE ($SPECIES, ${READS} reads in Feb): R1 files $r1, R2 files $r2, reference ${REF##*/}, $state"
        echo "  ${COMMAND[*]}"
    done
    exit 0
fi

count_command "$SLURM_ARRAY_TASK_ID"
if [ -e "$OUT/$LANE/outs/raw_feature_bc_matrix.h5" ]; then
    echo "$LANE: $OUT/$LANE/outs/raw_feature_bc_matrix.h5 exists, nothing to do"
    exit 0
fi
mkdir -p "$OUT/tmp/$LANE" "$CACHE"
export TMPDIR="$OUT/tmp/$LANE" XDG_CACHE_HOME="$CACHE" XDG_CONFIG_HOME="$CACHE"   # nothing goes to home
cd "$OUT"
echo "$LANE: started $(date) on $(hostname); expect ${READS} reads"
"${COMMAND[@]}"
echo "$LANE: finished $(date)"
grep -o '"Number of Reads","[^"]*"' "$OUT/$LANE/outs/web_summary.html" | head -1 || true
