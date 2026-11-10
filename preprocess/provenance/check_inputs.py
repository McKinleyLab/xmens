"""Checks the preprocessing inputs against records made when those inputs were created.

Reads samples.tsv, the lab's Samples sheet (CSV), the FASTQ folder and its md5 manifest, the February 2026 Cell Ranger
runs, both references and the mkref inputs of the mouse one, the transgene FASTA files, the files in sources/ and the
software named in config.sh. Writes a TSV with one row per check (per lane where it applies), each PASS, FAIL, PENDING
or INFO, prints a summary, and exits with status 1 if any check failed. provenance/check_inputs.sbatch gives the paths.
"""

import argparse
import gzip
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd


# Versions and reference releases stated in Methods.
CELLRANGER_VERSION = "cellranger-7.2.0"
CELLBENDER_VERSION = "0.2.2"
DOUBLETDETECTION_VERSION = "4.2"
HUMAN_REFERENCE_VERSION = "2024-A"  # 10x prebuilt GRCh38-2024-A
STOCK_MOUSE_REFERENCE_VERSION = "2020-A"  # 10x mm10-2020-A, the genome the transgenes were appended to

# Files mkref writes besides the STAR index in star/; anything else in a reference folder was put there later.
MKREF_FILES = {"reference.json", "fasta/genome.fa", "fasta/genome.fa.fai", "genes/genes.gtf.gz"}

READ_TYPES = ["I1", "I2", "R1", "R2"]
HASH_BLOCK_BYTES = 16 * 2**20
GTF_TAIL_BYTES = 2**16  # far more than the few transgene lines at the end of the GTF


def report(section: str, check: str, lane: str, result: str, observed: str, expected: str, evidence: str) -> dict:
    """One row of the report; result is PASS, FAIL, PENDING or INFO."""
    return {
        "section": section,
        "check": check,
        "lane": lane,
        "result": result,
        "observed": observed,
        "expected": expected,
        "evidence": evidence,
    }


def pass_fail(condition: bool) -> str:
    return "PASS" if condition else "FAIL"


def sha1_of(path: Path, n_bytes: int | None = None) -> str:
    """SHA-1 of a file, or of its first n_bytes; mkref and 10x record SHA-1 hashes in reference.json."""
    digest = hashlib.sha1()
    remaining = path.stat().st_size if n_bytes is None else n_bytes
    with path.open("rb") as handle:
        while remaining > 0:
            block = handle.read(min(HASH_BLOCK_BYTES, remaining))
            if not block:
                raise SystemExit(f"{path} ended {remaining:,} bytes early")
            digest.update(block)
            remaining -= len(block)
    return digest.hexdigest()


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Labels ----


def load_sheet(sheet_csv: Path) -> pd.DataFrame:
    """One row per lane of the Samples sheet, indexed by lane name as in samples.tsv (Mck 8a -> mck_8a)."""
    raw = pd.read_csv(sheet_csv, header=None, dtype=str)

    # Row 1 names the columns and row 2 gives units; lanes follow. The first column, the animal or donor ID, has no
    # name. An animal's ID, timepoint and age are written on its first lane only, so blanks below belong to it.
    names = raw.iloc[0].tolist()
    columns = {
        "animal": 0,
        "sample_timepoint": names.index("Sample timepoint"),
        "timepoint": names.index("Timepoint"),
        "age": names.index("Age @harvest"),
        "lane_id": names.index("Lane ID"),
        "gex_index": names.index("GEX Index"),
        "bcl": names.index("BCL object"),
    }
    sheet = raw.iloc[2:, list(columns.values())]
    sheet.columns = list(columns)
    sheet = sheet[sheet["lane_id"].str.startswith("Mck", na=False)].copy()

    animal_block = sheet["animal"].notna().cumsum()
    for column in ["animal", "sample_timepoint", "timepoint", "age"]:
        sheet[column] = sheet.groupby(animal_block)[column].ffill()

    # Human lanes have no Timepoint; their Sample timepoint is the cycle phase.
    sheet["condition"] = sheet["timepoint"].fillna(sheet["sample_timepoint"])
    sheet["lane"] = sheet["lane_id"].str.replace("Mck ", "mck_", regex=False)

    duplicated = sheet.loc[sheet["lane"].duplicated(), "lane"].tolist()
    if duplicated:
        raise SystemExit(f"{sheet_csv}: lanes listed twice: {duplicated}")
    return sheet.set_index("lane")


def plain(label: str) -> str:
    """Timepoint without underscores or case: samples.tsv writes 1_dpi_early where the sheet writes 1dpi_early."""
    return label.replace("_", "").lower()


def check_labels(samples: pd.DataFrame, sheet: pd.DataFrame, sheet_csv: Path) -> list[dict]:
    rows = []
    for sample in samples.itertuples(index=False):
        if sample.lane not in sheet.index:
            rows.append(report(
                section="labels",
                check="lane listed in the Samples sheet",
                lane=sample.lane,
                result="FAIL",
                observed="not listed",
                expected="listed",
                evidence=str(sheet_csv),
            ))
            continue

        entry = sheet.loc[sample.lane]
        rows.append(report(
            section="labels",
            check="GEX index (FASTQ prefix)",
            lane=sample.lane,
            result=pass_fail(sample.fastq_prefix == entry.gex_index),
            observed=sample.fastq_prefix,
            expected=entry.gex_index,
            evidence=f"Samples sheet, {entry.lane_id}",
        ))
        rows.append(report(
            section="labels",
            check="timepoint",
            lane=sample.lane,
            result=pass_fail(plain(sample.timepoint) == plain(entry.condition)),
            observed=sample.timepoint,
            expected=entry.condition,
            evidence=f"Samples sheet, {entry.lane_id}",
        ))
        rows.append(report(
            section="labels",
            check="animal or donor, age",
            lane=sample.lane,
            result="INFO",
            observed=f"{entry.animal}, {entry.age}",
            expected="-",
            evidence=f"Samples sheet, {entry.lane_id}",
        ))
    return rows


def lane_fastqs(fastq_root: Path, fastq_dir: str, prefix: str, read_type: str) -> list[Path]:
    # The underscore after the prefix keeps SI-TT-A1 from matching SI-TT-A10's files.
    return sorted((fastq_root / fastq_dir).glob(f"{prefix}_S*_L*_{read_type}_001.fastq.gz"))


def run_and_flowcell(fastq: Path) -> tuple[int, str]:
    """Run number and flowcell ID from the first read's header, @<instrument>:<run>:<flowcell>:<lane>:..."""
    with gzip.open(fastq, "rt") as handle:
        header = handle.readline()
    fields = header.lstrip("@").split(":")
    return int(fields[1]), fields[2]


def check_read_headers(samples: pd.DataFrame, sheet: pd.DataFrame, fastq_root: Path) -> list[dict]:
    """The sequencing run written into each lane's reads against the run (BCL folder) the sheet lists."""
    rows = []
    for sample in samples.itertuples(index=False):
        if sample.lane not in sheet.index:
            continue

        # BCL folders are named <date>_<instrument>_<run>_<flowcell>. NovaSeq X folders put the flowcell's side
        # (A or B) before its ID, which the read headers leave out; NextSeq folders have no side.
        bcl = sheet.loc[sample.lane, "bcl"]
        bcl_run, bcl_flowcell = bcl.split("_")[2:4]

        r1_files = lane_fastqs(fastq_root, sample.fastq_dir, sample.fastq_prefix, "R1")
        found = sorted({run_and_flowcell(path) for path in r1_files})
        matches = bool(found) and all(
            run == int(bcl_run) and bcl_flowcell.endswith(flowcell) for run, flowcell in found
        )
        rows.append(report(
            section="labels",
            check="sequencing run in the read headers",
            lane=sample.lane,
            result=pass_fail(matches),
            observed="; ".join(f"run {run}, flowcell {flowcell}" for run, flowcell in found) or "no R1 files",
            expected=bcl,
            evidence=f"first read of {len(r1_files)} R1 files",
        ))
    return rows


# FASTQs ----


def check_fastq_files(samples: pd.DataFrame, fastq_root: Path) -> list[dict]:
    rows = []
    for sample in samples.itertuples(index=False):
        n_files = {
            read_type: len(lane_fastqs(fastq_root, sample.fastq_dir, sample.fastq_prefix, read_type))
            for read_type in READ_TYPES
        }
        is_complete = n_files["R1"] > 0 and len(set(n_files.values())) == 1
        rows.append(report(
            section="FASTQs",
            check="complete set of files",
            lane=sample.lane,
            result=pass_fail(is_complete),
            observed=", ".join(f"{n} {read_type}" for read_type, n in n_files.items()),
            expected="the same number of I1, I2, R1 and R2 files",
            evidence=str(fastq_root / sample.fastq_dir),
        ))
    return rows


def check_fastq_manifest(fastq_root: Path, manifest: Path) -> list[dict]:
    """The md5 manifest lists exactly the FASTQs on disk; verifying the sums themselves is md5sum -c (hours)."""
    if not manifest.is_file() or manifest.stat().st_size == 0:
        return [report(
            section="FASTQs",
            check="md5 manifest lists every FASTQ",
            lane="-",
            result="PENDING",
            observed="no manifest yet",
            expected="one md5 line per FASTQ",
            evidence=str(manifest),
        )]

    # md5sum writes "<md5>  ./<folder>/<file>" when run from the FASTQ folder.
    listed = {line.split(maxsplit=1)[1].removeprefix("./") for line in manifest.read_text().splitlines() if line}
    on_disk = {str(path.relative_to(fastq_root)) for path in fastq_root.rglob("*.fastq.gz")}
    not_listed = sorted(on_disk - listed)
    not_on_disk = sorted(listed - on_disk)

    observed = f"{len(listed)} listed, {len(on_disk)} on disk"
    if not_listed or not_on_disk:
        observed += f"; not listed: {not_listed[:5]}; listed but missing: {not_on_disk[:5]}"
    return [report(
        section="FASTQs",
        check="md5 manifest lists every FASTQ",
        lane="-",
        result=pass_fail(not not_listed and not not_on_disk),
        observed=observed,
        expected="the same files",
        evidence=f"{manifest} (sha256 {sha256_of(manifest)[:16]})",
    )]


# Cell Ranger record ----


def invocation_field(invocation: str, pattern: str, path: Path) -> str:
    match = re.search(pattern, invocation)
    if match is None:
        raise SystemExit(f"{path}: no match for {pattern}")
    return match.group(1)


def check_cellranger_records(samples: pd.DataFrame, record_template: str, references: dict[str, Path]) -> list[dict]:
    """What each lane's February Cell Ranger run recorded about its inputs, settings and version."""
    rows = []
    for sample in samples.itertuples(index=False):
        pipestance = Path(record_template.format(lane=sample.lane))
        invocation_path = pipestance / "_invocation"
        invocation = invocation_path.read_text()

        sample_names = invocation_field(invocation, r'"sample_names": \[(.*?)\]', invocation_path).replace('"', "")
        read_path = invocation_field(invocation, r'"read_path": "([^"]*)"', invocation_path)
        reference_path = invocation_field(invocation, r'reference_path\s*=\s*"([^"]*)"', invocation_path)
        include_introns = invocation_field(invocation, r"include_introns\s*=\s*(\w+)", invocation_path)
        version = json.loads((pipestance / "_versions").read_text())["pipelines"]

        # The February run read decidualization/fastqs/<fastq_dir>; the FASTQs have since moved to the Dataverse
        # folder with the same <batch>/<index> layout. The references moved from /n/eddy_lab to
        # /n/lab_storage/eddy_lab, so only their folder names are compared.
        expected_read_path = f"/{sample.fastq_dir}"
        expected_reference = references[sample.species].name
        for check, is_expected, observed, expected in [
            ("sample index read", sample_names == sample.fastq_prefix, sample_names, sample.fastq_prefix),
            ("FASTQ folder read", read_path.endswith(expected_read_path), read_path, f"...{expected_read_path}"),
            (f"{sample.species} reference", Path(reference_path).name == expected_reference, reference_path,
             expected_reference),
            ("introns counted", include_introns == "true", include_introns, "true"),
            ("Cell Ranger version", version == CELLRANGER_VERSION, version, CELLRANGER_VERSION),
        ]:
            rows.append(report(
                section="Cell Ranger record",
                check=check,
                lane=sample.lane,
                result=pass_fail(is_expected),
                observed=observed,
                expected=expected,
                evidence=str(pipestance),
            ))
    return rows


# References ----


def fasta_name_and_length(path: Path) -> tuple[str, int]:
    """Name and sequence length of a single-record FASTA file."""
    lines = path.read_text().splitlines()
    headers = [line for line in lines if line.startswith(">")]
    if len(headers) != 1:
        raise SystemExit(f"{path}: expected one FASTA record, found {len(headers)}")
    name = headers[0][1:].split()[0]
    length = sum(len(line) for line in lines if not line.startswith(">"))
    return name, length


def as_fasta_record(path: Path) -> bytes:
    """A FASTA file's bytes as a record inside a larger FASTA, which ends in a newline even if the file doesn't."""
    record = path.read_bytes()
    return record if record.endswith(b"\n") else record + b"\n"


def last_lines(path: Path, n_lines: int) -> tuple[int, list[str]]:
    """Byte offset at which the last n_lines of a text file start, and those lines."""
    size = path.stat().st_size
    with path.open("rb") as handle:
        handle.seek(max(0, size - GTF_TAIL_BYTES))
        tail = handle.read()
    lines = tail.splitlines(keepends=True)[-n_lines:]
    offset = size - sum(len(line) for line in lines)
    return offset, [line.decode().rstrip("\n") for line in lines]


def describe_files_beyond_mkref(reference: Path) -> str:
    """Files mkref writes that are missing, and files it doesn't write that are present (star/ aside)."""
    files = {str(path.relative_to(reference)) for path in reference.rglob("*") if path.is_file()}
    added = sorted(name for name in files - MKREF_FILES if not name.startswith("star/"))
    missing = sorted(MKREF_FILES - files)
    return f"added: {', '.join(added) or 'none'}; missing: {', '.join(missing) or 'none'}"


def check_mouse_reference(
    built: Path,
    inputs: Path,
    transgene_fastas: list[Path],
    build_script: Path,
    build_script_record: Path,
) -> list[dict]:
    """The transgenic mouse reference against mkref's build record and 10x's record of stock mm10-2020-A."""
    build_record = json.loads((built / "reference.json").read_text())
    stock_record = json.loads((inputs / "reference.json").read_text())
    built_fasta = built / "fasta" / "genome.fa"
    input_fasta = inputs / "fasta" / "genome.fa"
    input_gtf = inputs / "genes" / "genes.gtf"
    stock_evidence = f"{inputs / 'reference.json'} (10x)"
    stock_name = f"10x mm10-{stock_record['version']}"

    # mkref records the SHA-1 of the genome.fa it wrote; its input should be the same file.
    built_fasta_sha1 = sha1_of(built_fasta)
    input_fasta_sha1 = sha1_of(input_fasta)

    # The input genome.fa should be 10x's stock genome.fa followed by the transgene records, byte for byte. The
    # transgene files end without a final newline; in genome.fa each record ends with one.
    transgene_records = b"".join(as_fasta_record(path) for path in transgene_fastas)
    n_stock_bytes = input_fasta.stat().st_size - len(transgene_records)
    stock_fasta_sha1 = sha1_of(input_fasta, n_bytes=n_stock_bytes)
    with input_fasta.open("rb") as handle:
        handle.seek(n_stock_bytes)
        appended = handle.read()

    # reference.json hashes the genes.gtf.gz that mkref wrote, which was later decompressed in place, so the
    # built genes.gtf is compared with mkref's input instead.
    built_gtf_sha1 = sha1_of(built / "genes" / "genes.gtf")
    input_gtf_sha1 = sha1_of(input_gtf)

    # The input GTF should be 10x's stock GTF followed by one line per transgene.
    transgenes = [fasta_name_and_length(path) for path in transgene_fastas]
    n_stock_gtf_bytes, transgene_lines = last_lines(input_gtf, len(transgenes))
    stock_gtf_sha1 = sha1_of(input_gtf, n_bytes=n_stock_gtf_bytes)

    checks = [
        ("genome.fa unchanged since mkref built it", built_fasta_sha1 == build_record["fasta_hash"],
         built_fasta_sha1, build_record["fasta_hash"], str(built / "reference.json")),
        ("mkref's input genome.fa is the built genome.fa", input_fasta_sha1 == built_fasta_sha1,
         input_fasta_sha1, built_fasta_sha1, str(input_fasta)),
        (f"genome.fa before the transgenes is {stock_name}",
         stock_fasta_sha1 == stock_record["fasta_hash"] and stock_record["version"] == STOCK_MOUSE_REFERENCE_VERSION,
         stock_fasta_sha1, stock_record["fasta_hash"], stock_evidence),
        ("genome.fa ends with exactly the transgene records", appended == transgene_records,
         f"last {len(appended):,} bytes", " + ".join(path.name for path in transgene_fastas), str(input_fasta)),
        ("genes.gtf identical to mkref's input", built_gtf_sha1 == input_gtf_sha1,
         built_gtf_sha1, input_gtf_sha1, str(input_gtf)),
        (f"genes.gtf before the transgenes is {stock_name}", stock_gtf_sha1 == stock_record["gtf_hash"],
         stock_gtf_sha1, stock_record["gtf_hash"], stock_evidence),
    ]

    # Each transgene is annotated as a single exon covering its whole sequence, on the + strand.
    for (name, length), line in zip(transgenes, transgene_lines, strict=True):
        fields = line.split("\t")
        observed = [fields[0], fields[2], fields[3], fields[4], fields[6]]
        expected = [name, "exon", "1", str(length), "+"]
        is_named = f'gene_id "{name}"' in fields[8] and f'transcript_id "{name}"' in fields[8]
        checks.append((f"{name} annotated as one exon spanning its sequence", observed == expected and is_named,
                       " ".join(observed), " ".join(expected), str(input_gtf)))

    # The build script in the reference folder, its copy in record/, and the mkref version it ran.
    script_text = build_script.read_text()
    mkref_line = next((line for line in script_text.splitlines() if " mkref " in line), "no mkref line")
    script_inputs = [f"--genome={built.name}", f"{inputs.name}/fasta/genome.fa", f"{inputs.name}/genes/genes.gtf"]
    mkref_version = build_record["mkref_version"].splitlines()[0]
    checks += [
        ("record/ copy of the build script is identical", sha256_of(build_script) == sha256_of(build_script_record),
         f"sha256 {sha256_of(build_script_record)[:16]}", f"sha256 {sha256_of(build_script)[:16]}",
         str(build_script_record)),
        ("build script runs mkref on these inputs", all(part in mkref_line for part in script_inputs),
         mkref_line.strip(), ", ".join(script_inputs), str(build_script)),
        ("mkref version", CELLRANGER_VERSION in mkref_version, mkref_version, CELLRANGER_VERSION,
         str(built / "reference.json")),
    ]

    rows = [
        report(
            section="mouse reference",
            check=check,
            lane="-",
            result=pass_fail(is_expected),
            observed=observed,
            expected=expected,
            evidence=evidence,
        )
        for check, is_expected, observed, expected, evidence in checks
    ]
    rows.append(report(
        section="mouse reference",
        check="files mkref didn't write",
        lane="-",
        result="INFO",
        observed=describe_files_beyond_mkref(built),
        expected="-",
        evidence=str(built),
    ))
    return rows


def check_human_reference(reference: Path) -> list[dict]:
    """The human reference against the hashes 10x recorded when it built GRCh38-2024-A."""
    record = json.loads((reference / "reference.json").read_text())
    evidence = f"{reference / 'reference.json'} (10x)"
    fasta_sha1 = sha1_of(reference / "fasta" / "genome.fa")
    gtf_sha1 = sha1_of(reference / "genes" / "genes.gtf.gz")

    checks = [
        ("genome.fa unchanged since 10x built it", fasta_sha1 == record["fasta_hash"], fasta_sha1,
         record["fasta_hash"]),
        ("genes.gtf.gz unchanged since 10x built it", gtf_sha1 == record["gtf_hash.gz"], gtf_sha1,
         record["gtf_hash.gz"]),
        ("release", record["version"] == HUMAN_REFERENCE_VERSION, record["version"], HUMAN_REFERENCE_VERSION),
    ]
    rows = [
        report(
            section="human reference",
            check=check,
            lane="-",
            result=pass_fail(is_expected),
            observed=observed,
            expected=expected,
            evidence=evidence,
        )
        for check, is_expected, observed, expected in checks
    ]
    rows.append(report(
        section="human reference",
        check="files mkref didn't write",
        lane="-",
        result="INFO",
        observed=describe_files_beyond_mkref(reference),
        expected="-",
        evidence=str(reference),
    ))
    return rows


# Software and sources ----


def installed_version(env: Path, package: str) -> str:
    """Version from the package's dist-info folder name; a set, because lib/python3.1 can link to lib/python3.10."""
    names = {path.name for path in env.glob(f"lib/python3*/site-packages/{package}-*.dist-info")}
    if len(names) != 1:
        return f"{len(names)} installs found"
    return names.pop().removeprefix(f"{package}-").removesuffix(".dist-info")


def check_software(cellranger: Path, cellbender_env: Path, qc_env: Path) -> list[dict]:
    version_line = subprocess.run(
        [str(cellranger), "--version"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    checks = [("Cell Ranger", version_line.endswith(CELLRANGER_VERSION), version_line, CELLRANGER_VERSION,
               str(cellranger))]

    for env, package, expected in [
        (cellbender_env, "cellbender", CELLBENDER_VERSION),
        (qc_env, "doubletdetection", DOUBLETDETECTION_VERSION),
    ]:
        version = installed_version(env, package)
        checks.append((package, version == expected, version, expected, str(env)))

    return [
        report(
            section="software",
            check=check,
            lane="-",
            result=pass_fail(is_expected),
            observed=observed,
            expected=expected,
            evidence=evidence,
        )
        for check, is_expected, observed, expected, evidence in checks
    ]


def check_sources(sources: Path) -> list[dict]:
    """Every file in sources/ is listed in SOURCES.tsv with the sha256 recorded when it was copied there."""
    listing = pd.read_csv(sources / "SOURCES.tsv", sep="\t")
    on_disk = sorted(
        str(path.relative_to(sources))
        for path in sources.rglob("*")
        if path.is_file() and path.name != "SOURCES.tsv"
    )

    rows = []
    for name in on_disk:
        recorded = listing.loc[listing["file"] == name, "sha256"].tolist()
        observed = sha256_of(sources / name)
        rows.append(report(
            section="sources",
            check="unchanged since copied",
            lane="-",
            result=pass_fail(recorded == [observed]),
            observed=f"{name}: sha256 {observed[:16]}",
            expected=f"sha256 {recorded[0][:16]}" if recorded else "listed in SOURCES.tsv",
            evidence=str(sources / "SOURCES.tsv"),
        ))
    for name in sorted(set(listing["file"]) - set(on_disk)):
        rows.append(report(
            section="sources",
            check="listed file present",
            lane="-",
            result="FAIL",
            observed=f"{name}: missing",
            expected="present",
            evidence=str(sources / "SOURCES.tsv"),
        ))
    return rows


# Run ----


def print_summary(report_table: pd.DataFrame, out: Path) -> None:
    counts = report_table.groupby(["section", "result"], sort=False).size().unstack(fill_value=0)
    print(f"\n{counts.to_string()}\n")

    for row in report_table[report_table["result"].isin(["FAIL", "PENDING"])].itertuples(index=False):
        print(f"{row.result} {row.section} | {row.check} | {row.lane}: observed {row.observed}; "
              f"expected {row.expected}")
    print(f"\nreport: {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--samples", type=Path, required=True, help="samples.tsv")
    parser.add_argument("--sheet", type=Path, required=True, help="the lab's Samples sheet as CSV")
    parser.add_argument("--fastqs", type=Path, required=True, help="the folder holding every lane's FASTQ folder")
    parser.add_argument("--fastq-md5", type=Path, required=True, help="md5sum output for every FASTQ under --fastqs")
    parser.add_argument("--cellranger-record", required=True,
                        help="a February Cell Ranger output folder, with {lane} in place of the lane")
    parser.add_argument("--ref-mouse", type=Path, required=True)
    parser.add_argument("--ref-mouse-inputs", type=Path, required=True,
                        help="the folder mkref built --ref-mouse from: stock mm10-2020-A plus the transgenes")
    parser.add_argument("--transgenes", type=Path, nargs="+", required=True,
                        help="the transgene FASTA files, in the order they were appended")
    parser.add_argument("--build-script", type=Path, required=True, help="the script that ran mkref")
    parser.add_argument("--build-script-record", type=Path, required=True, help="its copy in record/")
    parser.add_argument("--ref-human", type=Path, required=True)
    parser.add_argument("--cellranger", type=Path, required=True)
    parser.add_argument("--cellbender-env", type=Path, required=True)
    parser.add_argument("--qc-env", type=Path, required=True, help="the DoubletDetection environment")
    parser.add_argument("--sources", type=Path, required=True, help="provenance/sources, with SOURCES.tsv")
    parser.add_argument("--out", type=Path, required=True, help="report TSV; must not exist yet")
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f"{args.out} exists; not overwriting")

    samples = pd.read_csv(args.samples, sep="\t")
    sheet = load_sheet(args.sheet)
    references = {"mouse": args.ref_mouse, "human": args.ref_human}

    # The quick checks run first, so a label problem is in the log before the reference hashing starts.
    rows = []
    rows += check_labels(samples, sheet, args.sheet)
    rows += check_read_headers(samples, sheet, args.fastqs)
    rows += check_fastq_files(samples, args.fastqs)
    rows += check_fastq_manifest(args.fastqs, args.fastq_md5)
    rows += check_cellranger_records(samples, args.cellranger_record, references)
    rows += check_software(args.cellranger, args.cellbender_env, args.qc_env)
    rows += check_sources(args.sources)
    print(f"{len(rows)} label, FASTQ, record, software and source checks done; hashing the references", flush=True)

    rows += check_mouse_reference(
        args.ref_mouse,
        args.ref_mouse_inputs,
        args.transgenes,
        args.build_script,
        args.build_script_record,
    )
    rows += check_human_reference(args.ref_human)

    report_table = pd.DataFrame(rows)
    report_table.to_csv(args.out, sep="\t", index=False)
    print_summary(report_table, args.out)
    if (report_table["result"] == "FAIL").any():
        sys.exit(1)


if __name__ == "__main__":
    main()
