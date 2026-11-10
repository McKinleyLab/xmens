#!/usr/bin/env python3
"""Build deterministic human--mouse co-orthology graph components.

The output is an operational connected-component grouping derived from the
union of the two official Ensembl 116 species-specific homology dumps.  It is
not an official Ensembl HOG.  Every edge is joined to frozen stable gene-tree
content before graph construction.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence


RELEASE = "Ensembl 116"
HUMAN = "human"
MOUSE = "mouse"

ANNOTATION_COLUMNS = (
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

OUTPUT_COLUMNS = (
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

SUMMARY_COLUMNS = (
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

EDGE_COLUMNS = (
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

ATTRITION_COLUMNS = EDGE_COLUMNS + (
    "human_biotype",
    "mouse_biotype",
    "exclusion_reason",
)

EXPECTED_HOMOLOGY_TYPES = {
    "ortholog_one2one",
    "ortholog_one2many",
    "ortholog_many2many",
}


@dataclass(frozen=True)
class Gene:
    stable_id: str
    symbol: str
    biotype: str


@dataclass(frozen=True)
class Tree:
    homology_id: int
    root_id: int
    stable_id: str
    version: str
    member_type: str
    tree_type: str
    clusterset_id: str


@dataclass(frozen=True)
class Edge:
    human_gene: str
    mouse_gene: str
    human_protein: str
    mouse_protein: str
    homology_type: str
    homology_id: int
    human_identity: str
    mouse_identity: str
    goc_score: str
    wga_coverage: str
    is_high_confidence: str
    source_dump: str
    gene_tree_id: str

    def biological_key(self) -> tuple[object, ...]:
        return (
            self.human_gene,
            self.mouse_gene,
            self.human_protein,
            self.mouse_protein,
            self.homology_type,
            self.homology_id,
            self.human_identity,
            self.mouse_identity,
            self.goc_score,
            self.wga_coverage,
            self.is_high_confidence,
            self.gene_tree_id,
        )


Node = tuple[str, str]


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[Node, Node] = {}
        self.rank: dict[Node, int] = {}

    def add(self, item: Node) -> None:
        if item not in self.parent:
            self.parent[item] = item
            self.rank[item] = 0

    def find(self, item: Node) -> Node:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: Node, right: Node) -> None:
        self.add(left)
        self.add(right)
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    default_root = Path(__file__).resolve().parents[2]
    parser.add_argument("--repo-root", type=Path, default=default_root)
    parser.add_argument("--source-dir", type=Path, default=None,
                        help="frozen Ensembl 116 sources (default: <repo-root>/orthogroups/sources)")
    parser.add_argument("--out", type=Path, default=None,
                        help="component table to write (default: <repo-root>/orthogroups/...components.tsv)")
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_dict_rows(path: Path, expected_header: Sequence[str]) -> Iterator[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        observed = tuple(reader.fieldnames or ())
        if observed != tuple(expected_header):
            raise ValueError(f"unexpected header in {path}: {observed}")
        for row in reader:
            yield dict(row)


def load_annotation(path: Path, prefix: str) -> dict[str, Gene]:
    genes: dict[str, Gene] = {}
    for row in read_dict_rows(path, ANNOTATION_COLUMNS):
        stable_id = row["ensembl_gene_id"]
        if not stable_id.startswith(prefix):
            raise ValueError(f"unexpected gene ID in {path}: {stable_id}")
        if stable_id in genes:
            raise ValueError(f"duplicate gene ID in {path}: {stable_id}")
        genes[stable_id] = Gene(
            stable_id=stable_id,
            symbol=row["gene_symbol"],
            biotype=row["biotype"],
        )
    return genes


def load_trees(path: Path) -> dict[int, Tree]:
    trees: dict[int, Tree] = {}
    for row in read_dict_rows(path, TREE_COLUMNS):
        homology_id = int(row["homology_id"])
        tree = Tree(
            homology_id=homology_id,
            root_id=int(row["gene_tree_root_id"]),
            stable_id=row["gene_tree_id"],
            version=row["gene_tree_version"],
            member_type=row["member_type"],
            tree_type=row["tree_type"],
            clusterset_id=row["clusterset_id"],
        )
        if homology_id in trees:
            raise ValueError(f"duplicate gene-tree row for homology {homology_id}")
        if (
            not tree.stable_id
            or tree.member_type != "protein"
            or tree.tree_type not in {"tree", "supertree"}
            or tree.clusterset_id != "default"
        ):
            raise ValueError(f"unexpected stable gene-tree content: {tree}")
        trees[homology_id] = tree
    return trees


def normalize_edge(row: Mapping[str, str], source_dump: str, tree: Tree) -> Edge:
    left_species = row["species"]
    right_species = row["homology_species"]
    if {left_species, right_species} != {"homo_sapiens", "mus_musculus"}:
        raise ValueError(f"non-human--mouse row survived extraction: {row}")
    if left_species == "homo_sapiens":
        human_gene = row["gene_stable_id"]
        human_protein = row["protein_stable_id"]
        human_identity = row["identity"]
        mouse_gene = row["homology_gene_stable_id"]
        mouse_protein = row["homology_protein_stable_id"]
        mouse_identity = row["homology_identity"]
    else:
        mouse_gene = row["gene_stable_id"]
        mouse_protein = row["protein_stable_id"]
        mouse_identity = row["identity"]
        human_gene = row["homology_gene_stable_id"]
        human_protein = row["homology_protein_stable_id"]
        human_identity = row["homology_identity"]
    if not human_gene.startswith("ENSG") or not mouse_gene.startswith("ENSMUSG"):
        raise ValueError(f"invalid species gene IDs in homology {row['homology_id']}")
    homology_type = row["homology_type"]
    if homology_type not in EXPECTED_HOMOLOGY_TYPES:
        raise ValueError(f"unexpected homology type: {homology_type}")
    return Edge(
        human_gene=human_gene,
        mouse_gene=mouse_gene,
        human_protein=human_protein,
        mouse_protein=mouse_protein,
        homology_type=homology_type,
        homology_id=int(row["homology_id"]),
        human_identity=human_identity,
        mouse_identity=mouse_identity,
        goc_score=row["goc_score"],
        wga_coverage=row["wga_coverage"],
        is_high_confidence=row["is_high_confidence"],
        source_dump=source_dump,
        gene_tree_id=tree.stable_id,
    )


def swap_homology_orientation(row: Mapping[str, str]) -> dict[str, str]:
    """Return the same upstream homology row with its co-equal genes swapped."""
    swapped = dict(row)
    for left, right in (
        ("gene_stable_id", "homology_gene_stable_id"),
        ("protein_stable_id", "homology_protein_stable_id"),
        ("species", "homology_species"),
        ("identity", "homology_identity"),
    ):
        swapped[left], swapped[right] = row[right], row[left]
    return swapped


def load_edges(paths: Sequence[Path], trees: Mapping[int, Tree]) -> tuple[list[Edge], Counter[str]]:
    by_homology_id: dict[int, Edge] = {}
    source_counts: Counter[str] = Counter()
    source_names_by_id: defaultdict[int, set[str]] = defaultdict(set)
    for path in paths:
        for row in read_dict_rows(path, HOMOLOGY_COLUMNS):
            homology_id = int(row["homology_id"])
            if homology_id not in trees:
                raise ValueError(f"homology {homology_id} has no frozen stable tree row")
            edge = normalize_edge(row, path.name, trees[homology_id])
            swapped_edge = normalize_edge(
                swap_homology_orientation(row), path.name, trees[homology_id]
            )
            if edge.biological_key() != swapped_edge.biological_key():
                raise AssertionError(
                    f"orientation-dependent canonical edge for homology {homology_id}"
                )
            source_counts[path.name] += 1
            source_names_by_id[homology_id].add(path.name)
            previous = by_homology_id.get(homology_id)
            if previous is not None and previous.biological_key() != edge.biological_key():
                raise ValueError(f"conflicting duplicate homology ID {homology_id}")
            if previous is None:
                by_homology_id[homology_id] = edge
    if set(trees) != set(by_homology_id):
        extras = sorted(set(trees) - set(by_homology_id))
        raise ValueError(f"stable tree extract has homologies absent from union: {extras[:5]}")
    pair_to_id: dict[tuple[str, str], int] = {}
    for edge in by_homology_id.values():
        pair = (edge.human_gene, edge.mouse_gene)
        previous_id = pair_to_id.get(pair)
        if previous_id is not None and previous_id != edge.homology_id:
            raise ValueError(
                f"gene pair has multiple homology IDs: {pair} -> "
                f"{previous_id}, {edge.homology_id}"
            )
        pair_to_id[pair] = edge.homology_id
    return sorted(by_homology_id.values(), key=lambda edge: edge.homology_id), source_counts


def build_components(edges: Sequence[Edge]) -> list[tuple[tuple[str, ...], tuple[str, ...]]]:
    union_find = UnionFind()
    for edge in edges:
        union_find.union((HUMAN, edge.human_gene), (MOUSE, edge.mouse_gene))
    members: defaultdict[Node, list[Node]] = defaultdict(list)
    for node in union_find.parent:
        members[union_find.find(node)].append(node)
    result: list[tuple[tuple[str, ...], tuple[str, ...]]] = []
    for nodes in members.values():
        human = tuple(sorted(gene for species, gene in nodes if species == HUMAN))
        mouse = tuple(sorted(gene for species, gene in nodes if species == MOUSE))
        result.append((human, mouse))
    return sorted(result, key=lambda component: (component[0], component[1]))


def classify_component(
    human_count: int,
    mouse_count: int,
    homology_types: set[str],
    gene_tree_ids: set[str],
) -> str:
    if human_count == 0 or mouse_count == 0:
        return "incomplete"
    if not homology_types <= EXPECTED_HOMOLOGY_TYPES or len(gene_tree_ids) != 1:
        return "unusual"
    if human_count == 1 and mouse_count == 1:
        return "one_to_one" if homology_types == {"ortholog_one2one"} else "unusual"
    if human_count > 1 and mouse_count == 1:
        return (
            "human_expansion"
            if homology_types == {"ortholog_one2many"}
            else "unusual"
        )
    if human_count == 1 and mouse_count > 1:
        return (
            "mouse_expansion"
            if homology_types == {"ortholog_one2many"}
            else "unusual"
        )
    if human_count > 1 and mouse_count > 1:
        return (
            "both_expanded"
            if homology_types == {"ortholog_many2many"}
            else "unusual"
        )
    return "unusual"


def write_tsv_atomic(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> int:
    tmp_path = path.with_name(path.name + ".part")
    count = 0
    with tmp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
            count += 1
    os.replace(tmp_path, path)
    return count


def write_gzip_tsv_atomic(
    path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]
) -> int:
    tmp_path = path.with_name(path.name + ".part")
    count = 0
    with tmp_path.open("wb") as binary:
        with gzip.GzipFile(fileobj=binary, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                for row in rows:
                    writer.writerow(row)
                    count += 1
    os.replace(tmp_path, path)
    return count


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    source_dir = args.source_dir or args.repo_root / "orthogroups" / "sources"
    output_path = args.out or args.repo_root / "orthogroups" / "ensembl_human_mouse_components.tsv"

    human_annotation_path = source_dir / "ensembl_116_homo_sapiens_gene_annotation.tsv.gz"
    mouse_annotation_path = source_dir / "ensembl_116_mus_musculus_gene_annotation.tsv.gz"
    homology_paths = (
        source_dir / "ensembl_116_homo_sapiens_protein_default_human_mouse.tsv.gz",
        source_dir / "ensembl_116_mus_musculus_protein_default_human_mouse.tsv.gz",
    )
    tree_path = source_dir / "ensembl_116_human_mouse_gene_tree_content.tsv.gz"

    human_annotation = load_annotation(human_annotation_path, "ENSG")
    mouse_annotation = load_annotation(mouse_annotation_path, "ENSMUSG")
    trees = load_trees(tree_path)
    all_edges, source_counts = load_edges(homology_paths, trees)

    missing_human = sorted({edge.human_gene for edge in all_edges} - set(human_annotation))
    missing_mouse = sorted({edge.mouse_gene for edge in all_edges} - set(mouse_annotation))
    if missing_human or missing_mouse:
        raise ValueError(
            f"homology genes absent from frozen annotations: "
            f"human={missing_human[:5]}, mouse={missing_mouse[:5]}"
        )
    noncoding_edges = [
        edge
        for edge in all_edges
        if human_annotation[edge.human_gene].biotype != "protein_coding"
        or mouse_annotation[edge.mouse_gene].biotype != "protein_coding"
    ]
    edges = [edge for edge in all_edges if edge not in set(noncoding_edges)]

    # Ensembl documents arbitrary nonredundant placement between species files.
    # Completeness therefore means consuming both files and taking their union,
    # not expecting the same edge in each direction.
    if set(source_counts) != {path.name for path in homology_paths}:
        # Header-only inputs do not enter Counter; add them explicitly below.
        for path in homology_paths:
            source_counts.setdefault(path.name, 0)
    if set(source_counts) != {path.name for path in homology_paths}:
        raise AssertionError("both official species-specific dumps were not consumed")

    canonical_edges = {(edge.human_gene, edge.mouse_gene) for edge in edges}
    components = build_components(edges)
    node_to_component_index: dict[Node, int] = {}
    for index, (human_genes, mouse_genes) in enumerate(components, start=1):
        for gene in human_genes:
            if (HUMAN, gene) in node_to_component_index:
                raise AssertionError(f"human gene occurs in multiple components: {gene}")
            node_to_component_index[(HUMAN, gene)] = index
        for gene in mouse_genes:
            if (MOUSE, gene) in node_to_component_index:
                raise AssertionError(f"mouse gene occurs in multiple components: {gene}")
            node_to_component_index[(MOUSE, gene)] = index

    edges_by_component: defaultdict[int, list[Edge]] = defaultdict(list)
    for edge in edges:
        human_index = node_to_component_index[(HUMAN, edge.human_gene)]
        mouse_index = node_to_component_index[(MOUSE, edge.mouse_gene)]
        if human_index != mouse_index:
            raise AssertionError(f"edge crosses components: {edge.homology_id}")
        edges_by_component[human_index].append(edge)

    output_rows: list[tuple[object, ...]] = []
    summary_rows: list[tuple[object, ...]] = []
    component_classes: Counter[str] = Counter()
    multi_tree_components: list[str] = []
    incident_types: defaultdict[Node, set[str]] = defaultdict(set)
    incident_trees: defaultdict[Node, set[str]] = defaultdict(set)
    for edge in edges:
        for node in ((HUMAN, edge.human_gene), (MOUSE, edge.mouse_gene)):
            incident_types[node].add(edge.homology_type)
            incident_trees[node].add(edge.gene_tree_id)

    for index, (human_genes, mouse_genes) in enumerate(components, start=1):
        component_id = f"ENS116_HMCOMP_{index:06d}"
        component_edges = edges_by_component[index]
        homology_types = {edge.homology_type for edge in component_edges}
        gene_tree_ids = {edge.gene_tree_id for edge in component_edges}
        component_class = classify_component(
            len(human_genes), len(mouse_genes), homology_types, gene_tree_ids
        )
        component_classes[component_class] += 1
        if len(gene_tree_ids) != 1:
            multi_tree_components.append(component_id)
        all_coding = all(
            human_annotation[gene].biotype == "protein_coding" for gene in human_genes
        ) and all(
            mouse_annotation[gene].biotype == "protein_coding" for gene in mouse_genes
        )
        summary_rows.append(
            (
                component_id,
                component_class,
                len(human_genes),
                len(mouse_genes),
                len(component_edges),
                ";".join(sorted(homology_types)),
                ";".join(sorted(gene_tree_ids)),
                len(gene_tree_ids),
                "1" if all_coding else "0",
            )
        )
        for species, genes, annotation in (
            (HUMAN, human_genes, human_annotation),
            (MOUSE, mouse_genes, mouse_annotation),
        ):
            for gene_id in genes:
                gene = annotation[gene_id]
                node = (species, gene_id)
                output_rows.append(
                    (
                        component_id,
                        species,
                        gene_id,
                        gene.symbol,
                        gene.biotype,
                        ";".join(sorted(incident_types[node])),
                        ";".join(sorted(incident_trees[node])),
                        len(human_genes),
                        len(mouse_genes),
                        RELEASE,
                    )
                )

    if any(len(human) == 0 or len(mouse) == 0 for human, mouse in components):
        raise AssertionError("component lacking one of the focal species")
    if component_classes["incomplete"]:
        raise AssertionError("incomplete components were emitted")
    if component_classes["unusual"]:
        raise AssertionError(
            f"{component_classes['unusual']} components violate cardinality/tree invariants"
        )
    if any(
        len(human) == len(mouse) == 1
        and {edge.homology_type for edge in edges_by_component[index]}
        != {"ortholog_one2one"}
        for index, (human, mouse) in enumerate(components, start=1)
    ):
        raise AssertionError("a one-human/one-mouse component is not one-to-one")
    if any(
        edge.homology_type == "ortholog_one2many"
        and not (
            (len(components[node_to_component_index[(HUMAN, edge.human_gene)] - 1][0]) == 1)
            ^ (len(components[node_to_component_index[(HUMAN, edge.human_gene)] - 1][1]) == 1)
        )
        for edge in edges
    ):
        raise AssertionError("one-to-many edge has unexpected component copy numbers")
    if any(
        edge.homology_type == "ortholog_many2many"
        and not (
            len(components[node_to_component_index[(HUMAN, edge.human_gene)] - 1][0]) > 1
            and len(components[node_to_component_index[(HUMAN, edge.human_gene)] - 1][1]) > 1
        )
        for edge in edges
    ):
        raise AssertionError("many-to-many edge is not explicitly enumerated as both-expanded")

    output_count = write_tsv_atomic(output_path, OUTPUT_COLUMNS, output_rows)
    summary_path = source_dir / "ensembl_116_human_mouse_component_summary.tsv.gz"
    summary_count = write_gzip_tsv_atomic(summary_path, SUMMARY_COLUMNS, summary_rows)
    edge_path = source_dir / "ensembl_116_human_mouse_edges.tsv.gz"
    edge_rows = (
        (
            edge.human_gene,
            edge.mouse_gene,
            edge.human_protein,
            edge.mouse_protein,
            edge.homology_type,
            edge.homology_id,
            edge.gene_tree_id,
            edge.human_identity,
            edge.mouse_identity,
            edge.goc_score,
            edge.wga_coverage,
            edge.is_high_confidence,
            edge.source_dump,
        )
        for edge in edges
    )
    edge_count = write_gzip_tsv_atomic(edge_path, EDGE_COLUMNS, edge_rows)
    attrition_path = source_dir / "ensembl_116_non_protein_coding_homology_attrition.tsv.gz"
    attrition_rows = (
        (
            edge.human_gene,
            edge.mouse_gene,
            edge.human_protein,
            edge.mouse_protein,
            edge.homology_type,
            edge.homology_id,
            edge.gene_tree_id,
            edge.human_identity,
            edge.mouse_identity,
            edge.goc_score,
            edge.wga_coverage,
            edge.is_high_confidence,
            edge.source_dump,
            human_annotation[edge.human_gene].biotype,
            mouse_annotation[edge.mouse_gene].biotype,
            "one_or_both_GTF_gene_biotypes_not_protein_coding",
        )
        for edge in noncoding_edges
    )
    attrition_count = write_gzip_tsv_atomic(
        attrition_path, ATTRITION_COLUMNS, attrition_rows
    )

    excluded_biotype_pairs = Counter(
        (
            human_annotation[edge.human_gene].biotype,
            mouse_annotation[edge.mouse_gene].biotype,
        )
        for edge in noncoding_edges
    )

    validation = {
        "object_name": "ensembl_human_mouse_coortholog_component",
        "official_ensembl_hog": False,
        "release": RELEASE,
        "component_id_rule": (
            "lexicographically sort (sorted human gene tuple, sorted mouse gene tuple); "
            "assign ENS116_HMCOMP_%06d"
        ),
        "source_placement_note": (
            "Ensembl species-specific files contain arbitrary nonredundant subsets; "
            "both were consumed and unioned, so per-file mirror symmetry is not expected."
        ),
        "counts": {
            "human_annotation_genes": len(human_annotation),
            "mouse_annotation_genes": len(mouse_annotation),
            "source_rows": dict(sorted(source_counts.items())),
            "raw_unique_homology_ids": len(all_edges),
            "primary_protein_coding_homology_ids": len(edges),
            "excluded_non_protein_coding_homology_ids": len(noncoding_edges),
            "unique_canonical_gene_edges": len(canonical_edges),
            "stable_gene_tree_rows": len(trees),
            "components": len(components),
            "membership_rows": output_count,
            "component_summary_rows": summary_count,
            "edge_rows": edge_count,
            "non_protein_coding_attrition_rows": attrition_count,
            "excluded_biotype_pairs": {
                f"{human_biotype}|{mouse_biotype}": count
                for (human_biotype, mouse_biotype), count in sorted(
                    excluded_biotype_pairs.items()
                )
            },
            "human_member_genes": sum(len(human) for human, _ in components),
            "mouse_member_genes": sum(len(mouse) for _, mouse in components),
            "component_classes": dict(sorted(component_classes.items())),
            "homology_types": dict(sorted(Counter(e.homology_type for e in edges).items())),
            "multi_tree_components": len(multi_tree_components),
        },
        "assertions": {
            "both_species_specific_dumps_consumed": True,
            "union_homology_ids_unique_or_identical": True,
            "tree_join_complete": len(trees) == len(all_edges),
            "orientation_canonicalization_symmetric": True,
            "each_gene_in_one_component": len(node_to_component_index)
            == len({node for node in node_to_component_index}),
            "all_components_bipartite": all(h and m for h, m in components),
            "one_by_one_calls_are_one2one": True,
            "one2many_copy_numbers_valid": True,
            "many2many_explicitly_enumerated": True,
            "all_members_protein_coding": True,
            "one_stable_gene_tree_per_component": not multi_tree_components,
            "no_unusual_components": component_classes["unusual"] == 0,
        },
        "sha256": {
            output_path.name: sha256_file(output_path),
            summary_path.name: sha256_file(summary_path),
            edge_path.name: sha256_file(edge_path),
            attrition_path.name: sha256_file(attrition_path),
        },
        "multi_tree_component_ids": multi_tree_components,
    }
    validation_path = source_dir / "ensembl_116_component_validation.json"
    tmp_validation = validation_path.with_name(validation_path.name + ".part")
    with tmp_validation.open("w", encoding="utf-8") as handle:
        json.dump(validation, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp_validation, validation_path)

    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
