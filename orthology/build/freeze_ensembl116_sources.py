#!/usr/bin/env python3
"""Freeze compact Ensembl 116 inputs for the human--mouse orthogroup pilot.

This is the only Ensembl acquisition program in the pilot.  It downloads the
official release-116 GTF and species-specific protein_default homology files,
checks their upstream checksums, and writes deterministic compact extracts.
The two small SQL queries are read-only queries against Ensembl's public
release-116 MySQL databases.  Routine downstream analysis is fully offline.

PyMySQL is required only for acquisition (tested with PyMySQL 1.1.2):

    python -m pip install 'PyMySQL==1.1.2'
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sys
import tempfile
import urllib.request
from pathlib import Path
from typing import Iterable, Iterator, Sequence


RELEASE = "116"
MYSQL_HOST = "ensembldb.ensembl.org"
COMPARA_DB = "ensembl_compara_116"
HUMAN_CORE_DB = "homo_sapiens_core_116_38"
USER_AGENT = "human-ivd-orthogroup-pilot/1.0 (Ensembl release 116 freeze)"

GTF_COLUMNS = (
    "ensembl_gene_id",
    "gene_symbol",
    "biotype",
    "seqname",
    "start",
    "end",
    "strand",
)

HOMOLOGY_COLUMNS = (
    "gene_stable_id",
    "protein_stable_id",
    "species",
    "identity",
    "homology_type",
    "homology_gene_stable_id",
    "homology_protein_stable_id",
    "homology_species",
    "homology_identity",
    "dn",
    "ds",
    "goc_score",
    "wga_coverage",
    "is_high_confidence",
    "homology_id",
)

TREE_COLUMNS = (
    "homology_id",
    "gene_tree_root_id",
    "gene_tree_id",
    "gene_tree_version",
    "member_type",
    "tree_type",
    "clusterset_id",
)

MANIFEST_COLUMNS = (
    "source",
    "release_version",
    "retrieval_date",
    "url_endpoint",
    "query_parameters",
    "file_name",
    "sha256",
    "record_count",
    "taxonomic_level",
    "upstream_checksum",
)

ATTR_RE = re.compile(r'(\S+) "([^"]*)";')


SOURCES = {
    "human_gtf": {
        "url": (
            "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/"
            "Homo_sapiens.GRCh38.116.gtf.gz"
        ),
        "upstream_checksum": "bsd-sum:49151;blocks_1024:137815",
        "output": "ensembl_116_homo_sapiens_gene_annotation.tsv.gz",
        "species": "homo_sapiens",
    },
    "mouse_gtf": {
        "url": (
            "https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/"
            "Mus_musculus.GRCm39.116.gtf.gz"
        ),
        "upstream_checksum": "bsd-sum:42607;blocks_1024:105329",
        "output": "ensembl_116_mus_musculus_gene_annotation.tsv.gz",
        "species": "mus_musculus",
    },
    "human_homology": {
        "url": (
            "https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/"
            "homologies/homo_sapiens/"
            "Compara.116.protein_default.homologies.tsv.gz"
        ),
        "upstream_checksum": "md5:59857f48bbbdf6812999d58d7a24ccc4",
        "output": (
            "ensembl_116_homo_sapiens_protein_default_human_mouse.tsv.gz"
        ),
        "species": "homo_sapiens",
    },
    "mouse_homology": {
        "url": (
            "https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/"
            "homologies/mus_musculus/"
            "Compara.116.protein_default.homologies.tsv.gz"
        ),
        "upstream_checksum": "md5:8f9870f0f12ece5032f8e62117de9924",
        "output": (
            "ensembl_116_mus_musculus_protein_default_human_mouse.tsv.gz"
        ),
        "species": "mus_musculus",
    },
}

ENTREZ_SQL = """
SELECT DISTINCT g.stable_id AS ensembl_gene_id,
                x.dbprimary_acc AS entrez_id
FROM gene AS g
JOIN object_xref AS ox
  ON ox.ensembl_object_type = 'Gene'
 AND ox.ensembl_id = g.gene_id
JOIN xref AS x
  ON x.xref_id = ox.xref_id
JOIN external_db AS ed
  ON ed.external_db_id = x.external_db_id
WHERE ed.db_name = 'EntrezGene'
  AND g.stable_id LIKE 'ENSG%'
ORDER BY g.stable_id, x.dbprimary_acc
""".strip()

TREE_SQL_TEMPLATE = """
SELECT h.homology_id,
       h.gene_tree_root_id,
       r.stable_id AS gene_tree_id,
       r.version AS gene_tree_version,
       r.member_type,
       r.tree_type,
       r.clusterset_id
FROM homology AS h
LEFT JOIN gene_tree_root AS r
  ON r.root_id = h.gene_tree_root_id
WHERE h.homology_id IN ({placeholders})
ORDER BY h.homology_id
""".strip()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    default_root = Path(__file__).resolve().parents[2]
    parser.add_argument("--repo-root", type=Path, default=default_root)
    parser.add_argument(
        "--retrieval-date",
        default=dt.datetime.now(dt.timezone.utc).date().isoformat(),
        help="UTC retrieval date recorded in the metadata fragment",
    )
    parser.add_argument(
        "--keep-upstream",
        action="store_true",
        help="Keep downloaded upstream gzip files in the system temp directory",
    )
    return parser.parse_args()


def open_request(url: str):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(request, timeout=180)


def download(url: str, destination: Path) -> tuple[str, int]:
    sha256 = hashlib.sha256()
    size = 0
    with open_request(url) as response, destination.open("wb") as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            sha256.update(chunk)
            size += len(chunk)
    return sha256.hexdigest(), size


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bsd_sum_file(path: Path) -> tuple[int, int]:
    checksum = 0
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(chunk)
            for byte in chunk:
                checksum = ((checksum >> 1) | ((checksum & 1) << 15))
                checksum = (checksum + byte) & 0xFFFF
    return checksum, (size + 1023) // 1024


def verify_upstream(path: Path, declaration: str) -> None:
    fields = dict(item.split(":", 1) for item in declaration.split(";"))
    if "md5" in fields:
        observed = md5_file(path)
        if observed != fields["md5"]:
            raise RuntimeError(
                f"upstream MD5 mismatch for {path.name}: "
                f"expected {fields['md5']}, observed {observed}"
            )
        return
    observed_sum, observed_blocks = bsd_sum_file(path)
    expected_sum = int(fields["bsd-sum"])
    expected_blocks = int(fields["blocks_1024"])
    if (observed_sum, observed_blocks) != (expected_sum, expected_blocks):
        raise RuntimeError(
            f"upstream BSD checksum mismatch for {path.name}: expected "
            f"{expected_sum}/{expected_blocks}, observed "
            f"{observed_sum}/{observed_blocks}"
        )


def strip_version(stable_id: str) -> str:
    return re.sub(r"\.\d+$", "", stable_id)


def gzip_tsv_atomic(
    path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(path.name + ".part")
    count = 0
    with tmp_path.open("wb") as binary:
        with gzip.GzipFile(fileobj=binary, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                for row in rows:
                    writer.writerow(["" if value is None else value for value in row])
                    count += 1
    os.replace(tmp_path, path)
    return count


def extract_gtf(upstream: Path, output: Path) -> int:
    rows: list[tuple[str, str, str, str, int, int, str]] = []
    seen: set[str] = set()
    with gzip.open(upstream, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"malformed GTF row {line_number} in {upstream}")
            if fields[2] != "gene":
                continue
            attributes = dict(ATTR_RE.findall(fields[8]))
            stable_id = strip_version(attributes.get("gene_id", ""))
            biotype = attributes.get("gene_biotype", "")
            if not stable_id or not biotype:
                raise ValueError(
                    f"gene row {line_number} lacks gene_id/gene_biotype in {upstream}"
                )
            if stable_id in seen:
                raise ValueError(f"duplicate GTF gene ID {stable_id} in {upstream}")
            seen.add(stable_id)
            rows.append(
                (
                    stable_id,
                    attributes.get("gene_name", ""),
                    biotype,
                    fields[0],
                    int(fields[3]),
                    int(fields[4]),
                    fields[6],
                )
            )
    rows.sort(key=lambda row: row[0])
    return gzip_tsv_atomic(output, GTF_COLUMNS, rows)


def extract_human_mouse_homologies(upstream: Path, output: Path) -> int:
    rows: list[tuple[str, ...]] = []
    seen: set[tuple[str, ...]] = set()
    allowed_species = {"homo_sapiens", "mus_musculus"}
    with gzip.open(upstream, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = tuple(line.rstrip("\n").split("\t"))
            if fields == HOMOLOGY_COLUMNS:
                continue
            if len(fields) != len(HOMOLOGY_COLUMNS):
                raise ValueError(
                    f"homology row {line_number} has {len(fields)} fields, "
                    f"expected {len(HOMOLOGY_COLUMNS)} in {upstream}"
                )
            if {fields[2], fields[7]} != allowed_species:
                continue
            if fields in seen:
                raise ValueError(f"duplicate homology row in {upstream}: {fields}")
            seen.add(fields)
            rows.append(fields)
    rows.sort(key=lambda row: (int(row[14]), row[0], row[5], row[1], row[6]))
    return gzip_tsv_atomic(output, HOMOLOGY_COLUMNS, rows)


def require_pymysql():
    try:
        import pymysql  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "PyMySQL is required for acquisition only. Install exactly "
            "PyMySQL==1.1.2 and rerun."
        ) from exc
    return pymysql


def mysql_connection(database: str):
    pymysql = require_pymysql()
    return pymysql.connect(
        host=MYSQL_HOST,
        port=3306,
        user="anonymous",
        password="",
        database=database,
        charset="utf8mb4",
        connect_timeout=30,
        read_timeout=300,
        write_timeout=60,
    )


def freeze_entrez(output: Path) -> int:
    with mysql_connection(HUMAN_CORE_DB) as connection:
        with connection.cursor() as cursor:
            cursor.execute(ENTREZ_SQL)
            rows = [(strip_version(str(a)), str(b)) for a, b in cursor.fetchall()]
    if len(rows) != len(set(rows)):
        raise ValueError("duplicate Ensembl--Entrez pairs returned by Ensembl")
    if any(not gene.startswith("ENSG") or not entrez for gene, entrez in rows):
        raise ValueError("invalid Ensembl--Entrez row returned by Ensembl")
    return gzip_tsv_atomic(output, ("ensembl_gene_id", "entrez_id"), rows)


def homology_ids(paths: Sequence[Path]) -> list[int]:
    result: set[int] = set()
    for path in paths:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != HOMOLOGY_COLUMNS:
                raise ValueError(f"unexpected header in {path}")
            for row in reader:
                result.add(int(row["homology_id"]))
    return sorted(result)


def batched(values: Sequence[int], size: int) -> Iterator[Sequence[int]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def freeze_tree_content(homology_paths: Sequence[Path], output: Path) -> int:
    ids = homology_ids(homology_paths)
    rows: list[tuple[object, ...]] = []
    with mysql_connection(COMPARA_DB) as connection:
        with connection.cursor() as cursor:
            for batch in batched(ids, 750):
                placeholders = ",".join(["%s"] * len(batch))
                cursor.execute(
                    TREE_SQL_TEMPLATE.format(placeholders=placeholders), tuple(batch)
                )
                rows.extend(cursor.fetchall())
    returned_ids = [int(row[0]) for row in rows]
    if len(returned_ids) != len(set(returned_ids)):
        raise ValueError("duplicate homology IDs in stable gene-tree query")
    missing = sorted(set(ids) - set(returned_ids))
    if missing:
        raise ValueError(f"{len(missing)} homology IDs missing from Compara: {missing[:5]}")
    no_stable_tree = [row[0] for row in rows if row[2] in (None, "")]
    if no_stable_tree:
        raise ValueError(
            f"{len(no_stable_tree)} human--mouse homologies lack a stable tree ID"
        )
    rows.sort(key=lambda row: int(row[0]))
    return gzip_tsv_atomic(output, TREE_COLUMNS, rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest_fragment(path: Path, rows: list[dict[str, object]]) -> None:
    tmp_path = path.with_name(path.name + ".part")
    with tmp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(tmp_path, path)


def compact_sql(sql: str) -> str:
    return " ".join(sql.split())


def main() -> int:
    args = parse_args()
    source_dir = args.repo_root / "orthogroups" / "sources"
    source_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []

    temp_dir_obj = tempfile.TemporaryDirectory(prefix="ensembl116-freeze-")
    temp_dir = Path(temp_dir_obj.name)
    try:
        for key in ("human_gtf", "mouse_gtf", "human_homology", "mouse_homology"):
            spec = SOURCES[key]
            upstream = temp_dir / Path(str(spec["url"])).name
            print(f"Downloading {spec['url']}", file=sys.stderr, flush=True)
            download(str(spec["url"]), upstream)
            verify_upstream(upstream, str(spec["upstream_checksum"]))
            output = source_dir / str(spec["output"])
            if key.endswith("gtf"):
                count = extract_gtf(upstream, output)
                query = "GTF rows where feature == gene; retain named fields; sort by stable ID"
                level = "species annotation"
            else:
                count = extract_human_mouse_homologies(upstream, output)
                query = (
                    "retain all rows where unordered species pair == "
                    "{homo_sapiens,mus_musculus}; retain all 15 upstream fields; "
                    "sort by homology_id"
                )
                level = "human--mouse pairwise orthology"
            manifest_rows.append(
                {
                    "source": "Ensembl FTP",
                    "release_version": "Ensembl 116",
                    "retrieval_date": args.retrieval_date,
                    "url_endpoint": spec["url"],
                    "query_parameters": query,
                    "file_name": output.name,
                    "sha256": sha256_file(output),
                    "record_count": count,
                    "taxonomic_level": level,
                    "upstream_checksum": spec["upstream_checksum"],
                }
            )
            if args.keep_upstream:
                kept = temp_dir.with_name(temp_dir.name + "-kept")
                kept.mkdir(parents=True, exist_ok=True)
                shutil.copy2(upstream, kept / upstream.name)

        entrez_output = source_dir / "ensembl_116_homo_sapiens_entrez.tsv.gz"
        print("Querying release-116 human Entrez xrefs", file=sys.stderr, flush=True)
        entrez_count = freeze_entrez(entrez_output)
        manifest_rows.append(
            {
                "source": "Ensembl public MySQL",
                "release_version": "Ensembl 116; homo_sapiens_core_116_38",
                "retrieval_date": args.retrieval_date,
                "url_endpoint": f"mysql://{MYSQL_HOST}:3306/{HUMAN_CORE_DB}",
                "query_parameters": compact_sql(ENTREZ_SQL),
                "file_name": entrez_output.name,
                "sha256": sha256_file(entrez_output),
                "record_count": entrez_count,
                "taxonomic_level": "Homo sapiens (NCBI taxon 9606)",
                "upstream_checksum": "not applicable (versioned SQL query)",
            }
        )

        homology_paths = [
            source_dir / str(SOURCES["human_homology"]["output"]),
            source_dir / str(SOURCES["mouse_homology"]["output"]),
        ]
        tree_output = source_dir / "ensembl_116_human_mouse_gene_tree_content.tsv.gz"
        print("Querying stable gene-tree content", file=sys.stderr, flush=True)
        tree_count = freeze_tree_content(homology_paths, tree_output)
        manifest_rows.append(
            {
                "source": "Ensembl public MySQL",
                "release_version": "Ensembl 116; ensembl_compara_116",
                "retrieval_date": args.retrieval_date,
                "url_endpoint": f"mysql://{MYSQL_HOST}:3306/{COMPARA_DB}",
                "query_parameters": (
                    compact_sql(TREE_SQL_TEMPLATE).replace(
                        "{placeholders}",
                        "batched homology IDs from both frozen species extracts",
                    )
                ),
                "file_name": tree_output.name,
                "sha256": sha256_file(tree_output),
                "record_count": tree_count,
                "taxonomic_level": "human--mouse homology gene trees",
                "upstream_checksum": "not applicable (versioned SQL query)",
            }
        )

        manifest_path = source_dir / "ENSEMBL_SOURCE_MANIFEST_FRAGMENT.tsv"
        write_manifest_fragment(manifest_path, manifest_rows)
        summary = {
            "release": RELEASE,
            "retrieval_date": args.retrieval_date,
            "files": {row["file_name"]: row["record_count"] for row in manifest_rows},
        }
        print(json.dumps(summary, indent=2, sort_keys=True))
    finally:
        temp_dir_obj.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
