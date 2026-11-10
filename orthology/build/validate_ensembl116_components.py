#!/usr/bin/env python3
"""Independently validate the frozen Ensembl 116 component artifacts."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterator, Sequence


FINAL_HEADER = (
    "component_id",
    "species",
    "ensembl_gene_id",
    "gene_symbol",
    "biotype",
    "homology_type",
    "gene_tree_id",
    "human_copy_number",
    "mouse_copy_number",
    "source_release",
)
SUMMARY_HEADER = (
    "component_id",
    "component_class",
    "human_copy_number",
    "mouse_copy_number",
    "n_edges",
    "homology_types",
    "gene_tree_ids",
    "n_gene_tree_ids",
    "all_members_protein_coding",
)
EDGE_HEADER = (
    "human_ensembl_gene_id",
    "mouse_ensembl_gene_id",
    "human_protein_id",
    "mouse_protein_id",
    "homology_type",
    "homology_id",
    "gene_tree_id",
    "human_identity",
    "mouse_identity",
    "goc_score",
    "wga_coverage",
    "is_high_confidence",
    "source_dump",
)
ANNOTATION_HEADER = (
    "ensembl_gene_id",
    "gene_symbol",
    "biotype",
    "seqname",
    "start",
    "end",
    "strand",
)
RAW_HOMOLOGY_HEADER = (
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
TREE_HEADER = (
    "homology_id",
    "gene_tree_root_id",
    "gene_tree_id",
    "gene_tree_version",
    "member_type",
    "tree_type",
    "clusterset_id",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    return parser.parse_args()


def rows(path: Path, header: Sequence[str]) -> Iterator[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        observed = tuple(reader.fieldnames or ())
        if observed != tuple(header):
            raise AssertionError(f"unexpected header in {path}: {observed}")
        yield from reader


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expected_class(human: int, mouse: int, types: set[str], trees: set[str]) -> str:
    if not human or not mouse:
        return "incomplete"
    if len(trees) != 1:
        return "unusual"
    if human == mouse == 1:
        return "one_to_one" if types == {"ortholog_one2one"} else "unusual"
    if human > 1 and mouse == 1:
        return "human_expansion" if types == {"ortholog_one2many"} else "unusual"
    if human == 1 and mouse > 1:
        return "mouse_expansion" if types == {"ortholog_one2many"} else "unusual"
    if human > 1 and mouse > 1:
        return "both_expanded" if types == {"ortholog_many2many"} else "unusual"
    return "unusual"


def main() -> int:
    args = parse_args()
    source_dir = args.repo_root / "orthogroups" / "sources"
    final_path = args.repo_root / "orthogroups" / "ensembl_human_mouse_components.tsv"
    summary_path = source_dir / "ensembl_116_human_mouse_component_summary.tsv.gz"
    edge_path = source_dir / "ensembl_116_human_mouse_edges.tsv.gz"
    attrition_path = source_dir / "ensembl_116_non_protein_coding_homology_attrition.tsv.gz"
    tree_path = source_dir / "ensembl_116_human_mouse_gene_tree_content.tsv.gz"

    annotation: dict[tuple[str, str], tuple[str, str]] = {}
    for species, name in (
        ("human", "ensembl_116_homo_sapiens_gene_annotation.tsv.gz"),
        ("mouse", "ensembl_116_mus_musculus_gene_annotation.tsv.gz"),
    ):
        for row in rows(source_dir / name, ANNOTATION_HEADER):
            annotation[(species, row["ensembl_gene_id"])] = (
                row["gene_symbol"],
                row["biotype"],
            )

    members: defaultdict[str, dict[str, list[str]]] = defaultdict(
        lambda: {"human": [], "mouse": []}
    )
    gene_to_component: dict[tuple[str, str], str] = {}
    final_by_gene: dict[tuple[str, str], dict[str, str]] = {}
    final_count = 0
    for row in rows(final_path, FINAL_HEADER):
        final_count += 1
        species = row["species"]
        gene = row["ensembl_gene_id"]
        key = (species, gene)
        assert species in {"human", "mouse"}
        assert key not in gene_to_component
        assert key in annotation
        assert annotation[key] == (row["gene_symbol"], row["biotype"])
        assert row["biotype"] == "protein_coding"
        assert row["source_release"] == "Ensembl 116"
        gene_to_component[key] = row["component_id"]
        final_by_gene[key] = row
        members[row["component_id"]][species].append(gene)

    signatures = sorted(
        (
            tuple(sorted(by_species["human"])),
            tuple(sorted(by_species["mouse"])),
            component_id,
        )
        for component_id, by_species in members.items()
    )
    for index, (_, _, component_id) in enumerate(signatures, start=1):
        assert component_id == f"ENS116_HMCOMP_{index:06d}"

    summaries: dict[str, dict[str, str]] = {}
    classes: Counter[str] = Counter()
    for row in rows(summary_path, SUMMARY_HEADER):
        component_id = row["component_id"]
        assert component_id not in summaries
        assert component_id in members
        human_n = len(members[component_id]["human"])
        mouse_n = len(members[component_id]["mouse"])
        assert human_n == int(row["human_copy_number"])
        assert mouse_n == int(row["mouse_copy_number"])
        types = set(filter(None, row["homology_types"].split(";")))
        trees = set(filter(None, row["gene_tree_ids"].split(";")))
        assert int(row["n_gene_tree_ids"]) == len(trees) == 1
        assert row["component_class"] == expected_class(human_n, mouse_n, types, trees)
        assert row["component_class"] != "unusual"
        assert row["all_members_protein_coding"] == "1"
        for species in ("human", "mouse"):
            for gene in members[component_id][species]:
                final_row = final_by_gene[(species, gene)]
                assert int(final_row["human_copy_number"]) == human_n
                assert int(final_row["mouse_copy_number"]) == mouse_n
        summaries[component_id] = row
        classes[row["component_class"]] += 1
    assert set(summaries) == set(members)

    incident_types: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    incident_trees: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    primary_homology_ids: set[int] = set()
    primary_pairs: set[tuple[str, str]] = set()
    edge_count = 0
    edge_counts_by_component: Counter[str] = Counter()
    for row in rows(edge_path, EDGE_HEADER):
        edge_count += 1
        homology_id = int(row["homology_id"])
        pair = (row["human_ensembl_gene_id"], row["mouse_ensembl_gene_id"])
        assert homology_id not in primary_homology_ids
        assert pair not in primary_pairs
        primary_homology_ids.add(homology_id)
        primary_pairs.add(pair)
        human_key = ("human", pair[0])
        mouse_key = ("mouse", pair[1])
        assert human_key in gene_to_component and mouse_key in gene_to_component
        component_id = gene_to_component[human_key]
        assert component_id == gene_to_component[mouse_key]
        edge_counts_by_component[component_id] += 1
        for key in (human_key, mouse_key):
            incident_types[key].add(row["homology_type"])
            incident_trees[key].add(row["gene_tree_id"])
    assert set(gene_to_component) == set(incident_types) == set(incident_trees)
    for key, row in final_by_gene.items():
        assert row["homology_type"] == ";".join(sorted(incident_types[key]))
        assert row["gene_tree_id"] == ";".join(sorted(incident_trees[key]))
    for component_id, row in summaries.items():
        assert int(row["n_edges"]) == edge_counts_by_component[component_id]

    raw_ids: set[int] = set()
    raw_source_counts: dict[str, int] = {}
    for name in (
        "ensembl_116_homo_sapiens_protein_default_human_mouse.tsv.gz",
        "ensembl_116_mus_musculus_protein_default_human_mouse.tsv.gz",
    ):
        count = 0
        for row in rows(source_dir / name, RAW_HOMOLOGY_HEADER):
            homology_id = int(row["homology_id"])
            assert homology_id not in raw_ids
            raw_ids.add(homology_id)
            count += 1
        raw_source_counts[name] = count

    attrition_ids: set[int] = set()
    attrition_count = 0
    with gzip.open(attrition_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert tuple(reader.fieldnames or ())[: len(EDGE_HEADER)] == EDGE_HEADER
        for row in reader:
            attrition_count += 1
            homology_id = int(row["homology_id"])
            assert homology_id not in attrition_ids
            assert row["exclusion_reason"] == (
                "one_or_both_GTF_gene_biotypes_not_protein_coding"
            )
            assert row["human_biotype"] != "protein_coding" or row[
                "mouse_biotype"
            ] != "protein_coding"
            attrition_ids.add(homology_id)
    assert primary_homology_ids.isdisjoint(attrition_ids)
    assert primary_homology_ids | attrition_ids == raw_ids

    tree_ids: set[int] = set()
    for row in rows(tree_path, TREE_HEADER):
        homology_id = int(row["homology_id"])
        assert homology_id not in tree_ids
        assert row["gene_tree_id"]
        assert row["member_type"] == "protein"
        assert row["tree_type"] in {"tree", "supertree"}
        assert row["clusterset_id"] == "default"
        tree_ids.add(homology_id)
    assert tree_ids == raw_ids

    result = {
        "status": "PASS",
        "independent_assertions": {
            "exact_final_schema": True,
            "annotation_metadata_exact": True,
            "each_gene_in_exactly_one_component": True,
            "deterministic_component_ids_reproduced": True,
            "declared_copy_numbers_exact": True,
            "all_components_have_both_species": all(
                by_species["human"] and by_species["mouse"]
                for by_species in members.values()
            ),
            "component_classification_reproduced": True,
            "edge_endpoint_membership_exact": True,
            "incident_homology_types_exact": True,
            "incident_stable_tree_ids_exact": True,
            "raw_union_partitions_into_primary_plus_attrition": True,
            "stable_tree_join_complete": tree_ids == raw_ids,
        },
        "counts": {
            "components": len(members),
            "membership_rows": final_count,
            "primary_edges": edge_count,
            "non_protein_coding_attrition_edges": attrition_count,
            "raw_homologies": len(raw_ids),
            "stable_tree_rows": len(tree_ids),
            "classes": dict(sorted(classes.items())),
            "raw_source_rows": raw_source_counts,
        },
        "sha256": {final_path.name: sha256_file(final_path)},
    }
    output = source_dir / "ensembl_116_component_independent_validation.json"
    temporary = output.with_name(output.name + ".part")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, output)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
