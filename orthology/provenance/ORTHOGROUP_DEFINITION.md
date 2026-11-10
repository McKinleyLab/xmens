# Evolutionary object definition

## Primary object and verified taxon

The primary object is the set of extant human and mouse protein-coding genes
descending from **one ancestral gene in the human–mouse last common ancestral
lineage**. The operational taxonomic level is **Euarchontoglires (NCBI taxon
314146)**. NCBI Datasets and E-utilities independently identify 314146 as the
active superorder Euarchontoglires; OMA's own taxonomy places both Homo sapiens
and Mus musculus below the same node.

This level is suitable because duplications after that ancestral lineage yield
human or mouse co-orthologs inside one primary group, whereas duplications before
it produce separate primary groups. It therefore represents lineage-specific
copy expansion without merging older paralogs.

## A. Strict one-to-one ortholog pair

A strict one-to-one pair contains one human gene and one mouse gene related by
an orthology call whose cardinality is one-to-one. It is an edge-level relation,
not a synonym for a family. It supplies the conservative baseline.

## B. Human–mouse ancestral co-orthology class

A co-orthology class contains all human and mouse descendants of one
Euarchontoglires ancestral gene. It may be 1:1, 1:many, many:1, or many:many.
For Ensembl this pilot derives an operational
`ensembl_human_mouse_coortholog_component` as a connected component of the
official pairwise human–mouse orthology graph. These are **not official Ensembl
HOGs** and must be checked against reconciled gene trees because graph
transitivity can occasionally over- or under-group the intended ancestor.

## C. Taxon-specific hierarchical orthologous group (HOG)

An OMA HOG at Euarchontoglires is an explicitly inferred set of genes descending
from one ancestral gene at that named taxonomic level. Unlike the Ensembl graph
component, it is a resource-defined hierarchical object embedded in a species
tree and can have child HOGs after later duplications and a parent HOG before an
older duplication/speciation boundary.

## D. Broader parent HOG

A parent HOG represents a more ancient ancestral gene and can contain multiple
Euarchontoglires-level child HOGs created by duplications predating the
human–mouse ancestor. Parent identifiers are recorded where available, but
parent membership is not aggregated in this pilot. It is reserved for a later,
stricter analysis of redistribution among named child HOGs.

## E. Generic or root protein family

A root family is the broadest homologous/HOG collection under a deeply ancient
ancestor. It can join numerous old paralogs with different functions and is not
an appropriate primary conservation unit. Root identifiers are provenance only.

## Non-interchangeability

The five objects answer different evolutionary questions. A 1:1 pair tests
strict gene correspondence; a Euarchontoglires co-orthology class/HOG tests
descendants of one human–mouse ancestral gene; a parent HOG asks whether older
duplicate lineages redistribute activity; and a root family mostly asks about
broad homology. Substituting parent/root families for the primary object would
manufacture apparent rescue by pooling paralogs that were already distinct in
the human–mouse ancestor.

## Operational membership rules

- Primary computational key: unversioned Ensembl gene ID.
- Primary universe: protein-coding genes only. Noncoding DE rows stay visible in
  attrition and ordinary gene-level analyses.
- Ensembl components use only official release-116 human–mouse orthology edges;
  each component must contain both species and each gene may occur in only one
  component.
- OMA membership uses exactly the OMA HOG explicitly reported at
  Euarchontoglires. OMA protein IDs and every protein–Ensembl mapping are retained
  before any gene-level deduplication.
- Resource disagreements remain visible and receive manual tree review; neither
  source is silently preferred.

## Future mouse-inclusive states (specification only)

- `exact_one_to_one_conservation`: a strict 1:1 pair has eligible, same-sign
  responses under a predeclared effect/evidence rule.
- `coortholog_rescue`: the ancestral-group response has the same sign, but the
  dominant eligible human and mouse responders are not a strict 1:1 pair; the
  reconciled tree supports shared primary-group membership; weak/convergence-
  flagged genes do not drive the result; a size-matched null is insufficient;
  and Ensembl/OMA agree or a tree is manually adjudicated.
- `distributed_dosage_response`: multiple eligible paralogs collectively carry
  a same-direction response and no single gene dominates under a predeclared
  dominance threshold.
- `opposing_subfunctionalization`: eligible paralogs within at least one species
  show substantial opposing signed signal.
- `family_level_divergence`: comparable group-level responses disagree in sign
  or exceed a predeclared divergence threshold.
- `unresolved`: mapping, eligibility, resource disagreement, tree evidence, or
  power is insufficient for the other states.

With only two species these will be described as **candidate differential
paralog usage**, not definitive evolutionary paralog switches; additional
species are needed to polarize change.

## Older-parent-HOG plan (not executed)

For a later analysis, move only candidate cases to the nearest older parent HOG,
retain each Euarchontoglires child HOG identity, compare differential usage
among those child HOGs, and require stricter gene-tree and functional
adjudication. Never collapse indiscriminately to large root families.

