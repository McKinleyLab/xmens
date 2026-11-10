# Provenance of the preprocessing inputs

Every input the pipeline reads is checked against a record made when that input was created, by someone other than
this project wherever such a record exists: the lab's Samples sheet, the sequencer, 10x Genomics, mkref and Cell
Ranger. `check_inputs.sbatch` runs every check a script can do and writes one dated report. The items a script can't
check are listed under [Checked by a person](#checked-by-a-person).

The pipeline's outputs are checked separately: `compare_with_record.sbatch` compares a run with the run of record step
by step, and `notebooks/check_step_counts.ipynb` recounts the numbers in the shared objects.

## Run

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/preprocess
sbatch provenance/check_inputs.sbatch
```

It takes a few minutes, most of it hashing the two references (about 14 GB). The log,
`logs/check_inputs_<job ID>.out`, ends with a count of results per section and lists every FAIL and PENDING row. The
report is `provenance/reports/check_inputs_<date>_<time>.tsv`, one row per check, with the observed value, the expected
value and the file the expected value comes from. The job ends in FAILED if any check fails.

## What is checked, and against what

| Claim | Checked against | Recorded by |
|---|---|---|
| Each lane's index (FASTQ prefix) and timepoint in `samples.tsv` | `sources/samples_sheet_2026_11_06.csv` | the lab, when the samples were processed |
| Each lane's reads come from the sequencing run the sheet lists | run number and flowcell ID in the first read of every R1 file | the sequencer |
| Each lane has a complete set of FASTQs | I1, I2, R1 and R2 file counts | |
| The FASTQs are unchanged | `sources/fastq_md5_2026_11_06.txt` (once copied there; see below) | this project, 2026-11-06 |
| The February 2026 Cell Ranger run read each lane's index and FASTQ folder, used the species' reference, counted introns, and was version 7.2.0 | each lane's `_invocation` and `_versions` | Cell Ranger, February 2026 |
| The mouse genome.fa is unchanged since it was built | `fasta_hash` in `reference/transgenic_mus/reference.json` | mkref 7.2.0, 2024-03-01 |
| The mouse reference is stock mm10-2020-A with the two transgenes appended, nothing else | `fasta_hash` and `gtf_hash` in `reference/transgenic_mus_premkref/reference.json` | 10x Genomics, 2020 |
| The appended sequences are `sources/transgenes/*.fa`, each annotated as one exon spanning its sequence | the end of genome.fa and genes.gtf | |
| `mus.sh` built the mouse reference from those inputs, and `record/0_reference/mus.sh` is a copy of it | the script's mkref line | |
| The human reference is unchanged 10x GRCh38-2024-A | `fasta_hash` and `gtf_hash.gz` in its `reference.json` | 10x Genomics |
| Cell Ranger, CellBender and DoubletDetection versions are the ones in Methods | `cellranger --version`, installed package versions | |
| The files in `sources/` are unchanged since they were copied there | `sources/SOURCES.tsv` | this project |

## How the mouse reference was made

No script recorded the step before mkref; the checks above pin it down. `reference/transgenic_mus_premkref/` is 10x's
refdata-gex-mm10-2020-A with two files changed:

- `fasta/genome.fa`: the stock genome.fa, then `GqD_transgene.fa` and `GsD_Berdeaux.fa`, each followed by a newline
  (the files themselves end without one).
- `genes/genes.gtf`: the stock genes.gtf, then one line per transgene, `CHRM3_Jax026220` (3,592 bp) and
  `GsD_Berdeaux` (3,017 bp): a single exon on the + strand spanning the whole sequence, source `CAT`.

`mus.sh` then ran Cell Ranger 7.2.0 mkref on these two files on 2024-03-01, writing `reference/transgenic_mus/`.

## What these checks can't show

- **That today's FASTQs are the ones February's run read.** No checksums existed before 2026-11-06, and the February
  run read them from `decidualization/fastqs/`, which no longer exists. The link is reproduction: Cell Ranger outputs
  from today's FASTQs are byte-identical to February's in all 15 lanes (run 58258549;
  `compare_with_record.sbatch` repeats the comparison for any run).
- **That the STAR index files are intact.** Nothing recorded their hashes when they were built. The same
  reproduction covers them; rebuilding the reference with `mus.sh`'s command and comparing would too.
- **mkref's own hash of the mouse GTF.** `reference.json` hashes the `genes.gtf.gz` that mkref wrote; it has since
  been decompressed in place (the `genes/` folder changed on 2026-02-17), so the check compares the built
  `genes.gtf` with mkref's input instead. The report's INFO rows list such differences from mkref's output layout.

## Checked by a person

| Item | Status |
|---|---|
| The Samples sheet is the lab's record of which sample was loaded in each lane | Downloaded 2026-11-06 as CSV from the lab's "McKinley - SlideTags Samples" spreadsheet ("Samples" sheet). Who keeps it: not recorded yet |
| Where the two transgene sequences come from (plasmid map, accession or the lab that sent them) | Not recorded yet. The names point to JAX stock 026220 (GqD, `CHRM3_Jax026220`) and the Berdeaux lab (GsD) |
| Why these 15 of the sheet's sequenced lanes are used | The sheet notes cover some others (lanes clogged during GEM formation; poor spatial mapping noted on Mck 7a); not recorded for the rest |
| Methods matches the settings recorded here (introns, chemistry, versions, references) | Not checked yet |

## sources/

Copies of the documents the checks read, frozen on the date in their names; `SOURCES.tsv` records where each came
from and its sha256. Files here are never edited: a newer version goes in as a new dated file with its own row.

When the FASTQ checksum job (58355661) has finished, freeze its manifest here:

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/preprocess
MANIFEST=provenance/sources/fastq_md5_2026_11_06.txt
test ! -e "$MANIFEST" && cp -p /n/lab_storage/eddy_lab/users/nhilgert/tmp/2026_11_06_fastq_md5.txt "$MANIFEST"
printf '%s\t%s\t%s\t%s\n' fastq_md5_2026_11_06.txt \
    "md5sum of every FASTQ under xmens_final/dataverse/fastqs, job 58355661" 2026-11-06 \
    "$(sha256sum "$MANIFEST" | cut -d' ' -f1)" >> provenance/sources/SOURCES.tsv
```

`check_inputs.sbatch` then confirms the manifest lists exactly the FASTQs on disk. Recomputing every checksum reads
all 2.7 TB and takes hours, so it is a separate job; it prints only files whose checksum differs:

```bash
PIPELINE=/n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/preprocess
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/dataverse/fastqs
sbatch -p eddy -A eddy_lab -c 1 --mem 2G -t 12:00:00 -J fastq_md5_check -o "$PIPELINE/logs/%x_%j.out" \
    --wrap "md5sum -c --quiet $PIPELINE/provenance/sources/fastq_md5_2026_11_06.txt"
```

## When a check fails

The report row names the file the expected value came from. A FAIL means an input differs from its record. Find out
why before running the pipeline on it, and don't edit the record to make the check pass. PENDING means the record
doesn't exist yet (the FASTQ manifest until it is copied into `sources/`).
