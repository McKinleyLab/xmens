# Reproducing the February-regen QC'd objects

**Status (Nov 8, 2026):** checks 1–4 below pass; CellBender is running; QC, concat and positions are pending.

## The claim

The two QC'd objects sent to Cagri on Nov 6, 2026 can be rebuilt from the raw FASTQs with the scripts in this folder,
and every step gives the same result as the run that made them. Steps that are not bit-for-bit deterministic are named
in advance (acceptance criteria below), with the evidence for each.

## The target

The folder the message to Cagri pointed to, and the QC output its files were copied from (the copies only add the two
position sets). Pinned by sha256 on Nov 8 at 14:28, before any reproduction output past Cell Ranger existed
(`record/objects_sha256.txt`, Slurm job 58976537).

| File | Bytes | sha256 (first 16) | Nuclei |
|---|---|---|---|
| `2026_11_06_qc_regen_nick_andy_coords/qc_passed_mouse.h5ad` (sent) | 2,231,553,407 | `5f556b5e1603d3fc` | 69,307 |
| `2026_11_06_qc_regen_nick_andy_coords/qc_passed_human.h5ad` (sent) | 2,283,335,127 | `2849b88b10196a93` | 70,149 |
| `2026_11_06_qc_regen/qc_passed_mouse.h5ad` (QC output) | 2,222,286,119 | `fdae40f49d52ecb3` | 69,307 |
| `2026_11_06_qc_regen/qc_passed_human.h5ad` (QC output) | 2,273,794,314 | `a99d1b72d08f67f8` | 70,149 |

Both folders are in `/n/netscratch/eddy_lab/Everyone/nhilgert/`. Nuclei counts: `record/step_counts.tsv`, which also
gives the count at every step for every lane.

## How the target was made (the run of record)

Unedited copies of every script are in `record/`, each with its original path and sha256 in `record/SOURCES.tsv`
(all 34 checked against their originals on Nov 8). `record/LOGS_AND_OUTPUTS.md` says where every log is.

| Step | When | Script in `record/` | Output |
|---|---|---|---|
| 0 Mouse reference | 2024-03-01 | `0_reference/mus.sh` | `reference/transgenic_mus` |
| 1 Cell Ranger 7.2.0 | Feb 18 (mouse), Feb 28 (human) | `1_cellranger/` | `velocity/<lane>/cellranger/outs/` |
| 2 CellBender 0.2.2 | Jul 1 | `2_cellbender/01_cellbender.smk` | `xmens_inputs_regen/<sample>/` |
| 3 Ingest, QC, concat | Nov 6, 06:17–06:41 | `3_qc/` | `2026_11_06_qc_regen/` |
| 4 Positions | Nov 6, 07:24 | `4_positions/add_coords.py` | `2026_11_06_qc_regen_nick_andy_coords/` |

- Step 0 is one `cellranger mkref` 7.2.0 command. Step 1's exact command for each lane is in
  `1_cellranger/feb_cmdlines/`; its output folder is `decidualization/revisions/velocity/` on lab storage.
- Step 2 ran as Snakemake Slurm jobs 26739289–26739490: 8 samples on holygpu7c0530 and 7 on holygpu7c0522, each
  with one A40 GPU. Output: `decidualization/xmens_inputs_regen/` on lab storage.
- Step 3 used xmens_final commit 59e00dd, run as Slurm jobs 58262416 (15 samples) and 58262417 (concat), all on
  holygpu8a26505 (AMD Turin).

## How it is reproduced

`preprocess/` re-implements the run of record as five Slurm steps (see `README.md`), run into
`/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_xmens_preprocess/`:

| Step | Slurm job | Status |
|---|---|---|
| 1 Cell Ranger | 58329056 (Nov 6, 11:42–20:28) | done |
| 2 CellBender | 58978603 (Nov 8), all samples on holygpu7c0530 | running |
| 3 QC, 4 concat, 5 positions, comparison | queued behind step 2; QC restricted to AMD `genoa\|turin` nodes | pending |

## Checks

Every check is read-only and can be rerun by anyone. Each gives the commands, then what they showed on Nov 8.

### 0. Setup

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/preprocess
export OUT=/n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_06_xmens_preprocess
NS=/n/netscratch/eddy_lab/Everyone/nhilgert
SENT=$NS/2026_11_06_qc_regen_nick_andy_coords
QC_REC=$NS/2026_11_06_qc_regen
FEB=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/velocity
JULY=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/xmens_inputs_regen
REF=/n/lab_storage/eddy_lab/users/nhilgert/reference/transgenic_mus
```

### 1. The target files are the pinned ones

```bash
srun -p test -c 1 --mem=1G -t 00:20:00 sha256sum -c record/objects_sha256.txt
```

Expect `OK` for all four files. The file itself was written on Nov 8 by the same `sha256sum`, without `-c`.

### 2. The inputs are what the record says

```bash
jc=$(sbatch --parsable provenance/check_inputs.sbatch) && echo "check_inputs: $jc"
tail -30 logs/check_inputs_$jc.out                     # once it has finished, about 8 minutes
grep -w INFO provenance/reports/check_inputs_*.tsv | cut -c1-220
```

What it checks: the lab's Samples sheet against `samples.tsv`; each FASTQ's sequencing run and flowcell; the FASTQ
checksum list against the files on disk; the February Cell Ranger commands, versions and inputs; the frozen sources;
both references against the hashes Cell Ranger and 10x recorded.

**Nov 8** (Slurm job 58979047, 7 min 50 s; report `provenance/reports/check_inputs_2026_11_08_1446.tsv`):
174 checks, 157 PASS, 17 INFO, no FAIL, no PENDING.

| Section | PASS | INFO |
|---|---|---|
| labels | 45 | 15 |
| FASTQs | 16 | 0 |
| Cell Ranger record | 75 | 0 |
| software | 3 | 0 |
| sources | 4 | 0 |
| mouse reference | 11 | 1 |
| human reference | 3 | 1 |

INFO rows are facts for the reader, not checks:
- 15 labels rows give each lane's animal or donor and age from the Samples sheet: mice F355 (mck_8a, 8b), 1637 (mck_1),
  F565 (mck_12a, 12b), 5024 (mck_15), 2034 (mck_10a, 10b), F351 (mck_11a); donors D#6, 27 y (mck_6a–c, secretory) and
  D#16, 25 y (mck_5a–c, menstrual).
- Mouse reference: `genes/genes.gtf.gz` was decompressed in place to `genes/genes.gtf` (Feb 17, the day before the
  February mouse runs), and a `velocity.sh` script was saved in the folder. The GTF is identical to mkref's input
  (check 3 below); Cell Ranger does not read `velocity.sh`.
- Human reference: an Ensembl cDNA FASTA (`fasta/Homo_sapiens.GRCh38.cdna.all.fa`) was added to the 10x folder. It is
  not listed in the reference's `reference.json`, so `cellranger count` does not read it.

The FASTQ checksum list: `md5sum` of all 624 FASTQs (2.7 TB) under `xmens_final/dataverse/fastqs`, Slurm job 58355661
(Nov 6 13:38 – Nov 7 01:37), frozen as `provenance/sources/fastq_md5_2026_11_06.txt` (sha256 `f6802269ef0380b4`).
Recomputing all 624 checksums against it is a separate job of about 12 hours (see `provenance/README.md`); not yet run.

### 3. The mouse reference

```bash
for lane in $(cut -f2 samples.tsv | tail -n +2); do
    feb=$(grep -o -- '--transcriptome=[^ ]*' $FEB/$lane/cellranger/_cmdline | sed 's|.*/reference/||')
    new=$(grep -o -- '--transcriptome=[^ ]*' $OUT/cellranger/$lane/_cmdline | sed 's|.*/reference/||')
    echo "$lane  February: $feb  re-run: $new"
done
cat $REF/reference.json
grep mkref record/0_reference/mus.sh
tail -3 $REF/fasta/genome.fa.fai | cut -f1,2
tail -2 $REF/genes/genes.gtf | cut -f1-5,9
zcat $OUT/cellranger/mck_1/outs/filtered_feature_bc_matrix/features.tsv.gz | wc -l
zcat $OUT/cellranger/mck_1/outs/filtered_feature_bc_matrix/features.tsv.gz | tail -2
grep '^mouse reference' provenance/reports/check_inputs_2026_11_08_1446.tsv | cut -f2,4,5 | cut -c1-170
```

**Nov 8:**
- All 9 mouse lanes used `transgenic_mus` in February and in the re-run; all 6 human lanes used
  `human/refdata-gex-GRCh38-2024-A`. (February reached them through the old `/n/eddy_lab` mount; same folders.)
- `reference.json`: genome `transgenic_mus`, built by `cellranger mkref cellranger-7.2.0` from `genome.fa` and
  `genes.gtf`, `fasta_hash` `968a2ead03826a11ebea454a34e27536c72bfb9f`.
- `mus.sh` runs one `mkref` on `reference/transgenic_mus_premkref/fasta/genome.fa` and `.../genes/genes.gtf`.
- The genome ends with the last mm10 contig (`JH584295.1`, 1,976 bp) and two added contigs: `CHRM3_Jax026220`
  (GqD, 3,592 bp) and `GsD_Berdeaux` (GsD, 3,017 bp). Each is annotated as one exon spanning the contig, + strand.
- The re-run's mck_1 matrix lists 32,287 genes: 10x mm10-2020-A's 32,285 plus the two transgenes, listed last.
- The 11 PASS checks:

| Check | Value |
|---|---|
| genome.fa unchanged since mkref built it | `968a2ead…` |
| mkref's input genome.fa is the built genome.fa | `968a2ead…` |
| genome.fa before the transgenes is 10x mm10-2020-A | `c9c31fef…` |
| genome.fa ends with exactly the transgene records | last 6,642 bytes |
| genes.gtf identical to mkref's input | `0c7b0e5a…` |
| genes.gtf before the transgenes is 10x mm10-2020-A | `a1e5c293…` |
| CHRM3_Jax026220 annotated as one exon spanning its sequence | exon 1–3592, + |
| GsD_Berdeaux annotated as one exon spanning its sequence | exon 1–3017, + |
| record/ copy of the build script is identical | sha256 `dda5c508…` |
| build script runs mkref on these inputs | the `mkref` command above |
| mkref version | cellranger-7.2.0 |

### 4. Cell Ranger: a fresh run, byte-identical to February

```bash
cat $OUT/cellranger/mck_1/_versions; echo
tr ' ' '\n' < $OUT/cellranger/mck_1/_cmdline | grep -E -- '--(id|fastqs|transcriptome|sample)'
tr ' ' '\n' < $FEB/mck_1/cellranger/_cmdline | grep -E -- '--(id|fastqs|transcriptome|sample)'
ls -l --time-style=long-iso $OUT/cellranger/mck_1/outs/raw_feature_bc_matrix.h5 \
    $FEB/mck_1/cellranger/outs/raw_feature_bc_matrix.h5
for lane in $(cut -f2 samples.tsv | tail -n +2); do
    printf '%-8s' $lane
    for f in raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5 metrics_summary.csv; do
        if cmp -s $OUT/cellranger/$lane/outs/$f $FEB/$lane/cellranger/outs/$f; then printf ' same'
        else printf ' DIFFERENT(%s)' $f; fi
    done; echo
done
```

**Nov 8:**
- Versions: `cellranger-7.2.0`, Martian `v4.0.11`.
- Same lane, sample index (`SI-TT-A5`) and reference. Different FASTQ location: February read
  `decidualization/lanes/mck_1/gex` (now an empty folder); the re-run read
  `xmens_final/dataverse/fastqs/mck_1/SI-TT-A5`, the files the checksum list pins.
- A new file, not a copy: mck_1's raw matrix is 48,847,180 bytes in both, dated Feb 18 15:55 and Nov 6 12:51.
- **All 15 lanes: raw matrix, filtered matrix and `metrics_summary.csv` byte-identical.** The metrics include each
  lane's read count, mapping rates and sequencing saturation, so identical files also show the same reads went in.

### 5. CellBender against the July run

```bash
# Barcode setup, decided from the raw matrix before training (one sample)
grep -h -E 'probable cell barcodes|surely-empty' $JULY/uninduced_mck8a/uninduced_mck8a_cellbender.log \
    $OUT/logs/2_cellbender_uninduced_mck8a.log
# Cells called, every finished sample
for d in $OUT/cellbender/*/; do
    s=$(basename $d); f=${s}_cellbender_cell_barcodes.csv
    [ -e $d/$f ] || continue
    cmp -s $d/$f $JULY/$s/$f && echo "$s: same cells as July" || echo "$s: DIFFERENT cells from July"
done
# Counts per cell: compare_with_record (end of the chain, or any time with: sbatch compare_with_record.sbatch)
```

**Nov 8:** barcode setup identical in `uninduced_mck8a`: 13,007 probable cells plus 36,993 further barcodes (the
50,000 of `--total-droplets-included`), 45,724 empty droplets, largest surely-empty droplet 215 UMIs.
**Cells called and counts: pending.**

### 6. QC numbers against the record's logs

```bash
for s in $(cut -f1 samples.tsv | tail -n +2); do
    pattern='loaded [0-9]+ cells|after threshold QC: [0-9]+|flagged [0-9]+ doublets|final: [0-9]+ cells'
    a=$(grep -h -o -E "$pattern" $QC_REC/logs/01_qc_$s.log)
    b=$(grep -h -o -E "$pattern" $OUT/logs/3_qc_$s.log)
    [ "$a" = "$b" ] && echo "$s: same QC numbers" || echo "$s: DIFFERENT QC numbers"
done
```

The pattern picks out four numbers per sample (cells loaded, passing thresholds, doublets flagged, final). Tested on
Nov 8: every record log yields exactly four, and for mck_1 they are 4,106, 3,567, 70 and 3,497, as in
`record/step_counts.tsv`. On the one-sample pipeline test of Nov 6 (`2026_11_06_preprocess_test/`) it finds 68
doublets and 3,499 final nuclei, the Intel result in the table further down, so the check would report that sample
as different.

**Pending.**

### 7–8. Final objects and positions

```bash
tail -40 logs/compare_with_record_*.out        # the summary
ls $OUT/compare_with_record_*.tsv              # one row per sample, every step
```

`compare_with_record` checks, per sample: CellBender cells and identical counts against July; the QC numbers and the
whole QC log against the record's; nuclei kept and identical counts against `2026_11_06_qc_regen/qc_passed_*.h5ad`;
`_nick` and `_andy` positions on shared nuclei against the files sent to Cagri. **Pending.**

## Acceptance criteria (proposed; fixed before the results)

| Step | Reproduced if | Nov 8 |
|---|---|---|
| Cell Ranger | raw and filtered matrices and metrics byte-identical, 15 of 15 lanes | met |
| CellBender | same called cells and identical counts for every cell, 15 of 15 samples | pending |
| QC | identical loaded, threshold, doublet and final numbers, 15 of 15 samples | pending |
| Concat | same nuclei in each `qc_passed` object, identical counts for every nucleus | pending |
| Positions | identical `_nick` and `_andy` positions for every nucleus | pending |

Any shortfall is reported against these criteria, with its size and cause; the criteria are not changed afterward.

## Earlier evidence on the two non-deterministic steps

**CellBender** trains on the GPU. Two runs on mck_1 (June 30 and July 1, 2026, both on an A40 on holygpu7c0530 with
the same environment) agreed in every field of the output: counts, latent variables and training curve. The files
differ only in HDF5 metadata bytes, so contents are compared, not checksums.

**DoubletDetection** (in QC) depends on the CPU and thread count, not on chance: on mck_1 (Nov 6),

| Run | Node | CPU | Threads | Doublets | Kept nuclei |
|---|---|---|---|---|---|
| record | holygpu8a26505 | AMD Turin | 8 | 70 | set A |
| test | holy8a28511 | AMD Genoa | 8 | 70 | identical to A |
| test | holy8a27111 | Intel Ice Lake | 8 | 68 | set B |
| test | holy8a27312 | Intel Ice Lake | 8 | 68 | identical to B |
| two tests | holy8a27211 | Intel Ice Lake | 1 | 74 | set C, identical in both |

Thresholds were identical in every run. Caveats: one sample, two nodes per setting, and CPU type is confounded with core
count (Intel nodes 64 cores, AMD 192). Hence step 3 runs on AMD `genoa|turin` nodes.

**QC scripts on the 2025 inputs** (Nov 5, `2026_11_05_qc_verify_scripts/`): ingest, threshold QC and the QC objects on
shared nuclei identical in 15 of 15 samples; doublet calls differed in all 15 (26,311 vs 26,444 flagged).

## Differences from the run of record

Each difference, and why it does not change the result:

- **Paths:** `/n/eddy_lab/...` became `/n/lab_storage/eddy_lab/...`. The same files under a new mount.
- **FASTQs:** read from `xmens_final/dataverse/fastqs`, because February's `decidualization/lanes/<lane>/gex` folders
  are now empty. The checksum list pins these files, and Cell Ranger's output is byte-identical (check 4).
- **Cell Ranger flags:** `--no-bam` and `--disable-ui` are added, and step 1 avoids the nodes the February script
  excluded. Matrices and metrics are byte-identical (check 4).
- **Reference folders:** the mouse folder holds a decompressed `genes.gtf` (Feb 17, before the February runs) and
  `velocity.sh`; the human folder holds an extra cDNA FASTA. The GTF is identical to mkref's input, the other two are
  not read, and the output is byte-identical (checks 3 and 4).
- **CellBender launcher:** run as `$CELLBENDER_ENV/bin/python $CELLBENDER_ENV/bin/cellbender`. The launcher's first line
  names the old `/n/eddy_lab` interpreter path, so the first attempt (job 58976648) failed with "bad interpreter".
  This is the same four-line launcher on the same interpreter.
- **CellBender hardware:** all 15 samples run on holygpu7c0530. In July, 8 ran there and 7 on holygpu7c0522, which has
  been down since Nov 3; both nodes have 8 A40 GPUs and Intel Cascade Lake CPUs. Tested by check 5.
- **QC scripts:** Snakemake inputs became command-line arguments, and `qc.py` stops if DoubletDetection is missing or
  fails, instead of warning. The computations are unchanged; tested by checks 6–8.

## Known gaps

- **Mouse reference:** no script for adding the transgene contigs to the FASTA and GTF, and no mkref log. Bounded by
  check 3: everything before the transgenes is 10x's stock files, and the last 6,642 bytes are the transgene records.
- **QC config:** the record's QC manifests say `config/config.yaml` had uncommitted changes relative to commit 59e00dd.
  The file as run is `record/config/config.yaml` (sha256 in `record/SOURCES.tsv`); the diff is not yet written down.
- **CellBender logs:** the July Slurm logs are gone; each sample's own log keeps its command line.
- **Andy's positions:** no scripts or logs exist; they are an external input and need checksums.
- **On netscratch only:** the target objects, the record QC run's logs, manifests and specs, its Slurm logs, and this
  reproduction's run folder. Netscratch deletes files after 90 days.
- **Not in `record/`:** the three comparison scripts that `run_ingest_qc.sbatch` calls.
- **Environments:** `cellbender_022`, `cuda11.8_env` and `rsc_uv` are not yet exported into the repo.

## Still to do

1. Finish steps 2–5; record checks 5–8 here.
2. Copy the target, the record QC run's logs and this run's results to lab storage; point `compare_with_record` at
   the copies.
3. Write down the `config.yaml` diff; add checksums for both position sets; copy the three comparison scripts.
4. Export the three environments.
5. Run the 12-hour FASTQ re-hash once.
6. Commit and tag; cite the tag here.
