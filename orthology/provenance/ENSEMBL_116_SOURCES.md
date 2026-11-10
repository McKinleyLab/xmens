# Ensembl 116 frozen-source contract

These files pin every Ensembl input used by the human--mouse operational
co-orthology layer to **Ensembl release 116**. Acquisition is the only live
step. Component construction reads only the frozen local files.

## Acquisition and build

The acquisition script is
`orthogroup_work/scripts/freeze_ensembl116_sources.py`. It was run with
PyMySQL 1.1.2 and retrieval date `2026-08-07`. It verifies the official FTP
checksum before extracting either GTF or homology data. The offline builder is
`orthogroup_work/scripts/build_ensembl116_components.py`.

The external-input metadata fragment is
`ENSEMBL_SOURCE_MANIFEST_FRAGMENT.tsv`. It is intentionally not named
`SOURCE_MANIFEST.tsv`, because the coordinating analysis combines this fragment
with independently acquired OMA and taxonomy sources.

## Frozen external extracts

`ensembl_116_homo_sapiens_gene_annotation.tsv.gz` and
`ensembl_116_mus_musculus_gene_annotation.tsv.gz` contain one official GTF
`gene` feature per row, sorted by unversioned stable gene ID. Their exact schema
is:

```text
ensembl_gene_id  gene_symbol  biotype  seqname  start  end  strand
```

The annotation universe is Ensembl's standard release-116 GTF for GRCh38 or
GRCm39. It is not the separate `chr_patch_hapl_scaff` GTF.

`ensembl_116_homo_sapiens_entrez.tsv.gz` contains every distinct gene-level
`EntrezGene` xref attached to an ENSG record in the official
`homo_sapiens_core_116_38` database:

```text
ensembl_gene_id  entrez_id
```

This lossless core-xref extract includes 5,074 pairs representing 3,787 ENSG
IDs absent from the frozen standard GTF (core records can extend beyond that
GTF universe). Identifier
harmonization must intersect it with the frozen GTF and report those excluded
pairs; it must not silently treat every core xref as a current pilot target.

The two species-specific `protein_default` extracts preserve all 15 upstream
homology fields for the unordered species pair `{homo_sapiens,mus_musculus}`:

```text
gene_stable_id  protein_stable_id  species  identity  homology_type
homology_gene_stable_id  homology_protein_stable_id  homology_species
homology_identity  dn  ds  goc_score  wga_coverage  is_high_confidence
homology_id
```

Ensembl explicitly documents that its species-specific files contain arbitrary
nonredundant subsets. Both species files are therefore required and their union,
not an intersection, is the complete pairwise input. In release 116 the human
extract contributes zero human--mouse rows and the mouse extract contributes all
23,764. This source asymmetry is expected and retained explicitly.

`ensembl_116_human_mouse_gene_tree_content.tsv.gz` joins every one of those
23,764 homology IDs to the official stable tree root in
`ensembl_compara_116`:

```text
homology_id  gene_tree_root_id  gene_tree_id  gene_tree_version
member_type  tree_type  clusterset_id
```

The file contains 23,286 ordinary `tree` rows and 478 `supertree` rows, all in
the `protein/default` collection. It is sufficient to place every operational
edge in a stable gene-tree object. It does **not** freeze full topology or
reconciled internal duplication-node attributes; targeted tree inspection is
still required for manual duplication-timing adjudication.

## Derived audit files

- `ensembl_116_human_mouse_edges.tsv.gz` is the canonical, orientation-neutral,
  protein-coding edge table used by the graph.
- `ensembl_116_non_protein_coding_homology_attrition.tsv.gz` retains all 1,072
  excluded edges and both endpoint biotypes. Most are IG/TR gene biotypes.
- `ensembl_116_human_mouse_component_summary.tsv.gz` contains one row per
  component with its copy-number class, edge count, homology types, and stable
  tree IDs.
- `ensembl_116_component_validation.json` records counts, SHA256 checksums, and
  every automated invariant.

The final membership table has exactly the requested schema:

```text
component_id  species  ensembl_gene_id  gene_symbol  biotype
homology_type  gene_tree_id  human_copy_number  mouse_copy_number
source_release
```

`homology_type` and `gene_tree_id` are sorted semicolon-delimited sets of the
gene's incident edges. Component identifiers are deterministic: components are
sorted lexicographically by their sorted human and mouse member tuples and then
numbered `ENS116_HMCOMP_000001`, ... . These are operational
`ensembl_human_mouse_coortholog_component` IDs, **not official Ensembl HOGs**.
