# orthology: the human-mouse ortholog table (Ensembl 116)

The table of human-mouse orthologs behind the cross-species gene comparisons, with the code and the frozen sources
to rebuild it.

| file | contents |
|---|---|
| `ensembl_human_mouse_components.tsv` | 36,334 protein-coding genes with a human-mouse homolog, in 17,030 components |
| `ortholog_pairs_one2one_gene_ids.csv` | the 16,314 one-to-one pairs by Ensembl gene ID, with symbols |
| `ortholog_pairs_one2one.csv` | the 16,275 of those pairs whose two genes each have a symbol of their own |
| `orthology.py` | a module to query the table (`one_to_one`, `human_to_mouse`, `mouse_to_human`, `homology`) |
| `check_orthology.py` | checks this folder against its checksums and the other copies on the lab's storage |
| `build/` | the builder and the frozen Ensembl 116 sources; `bash build/rebuild.sh` rebuilds the table |
| `provenance/` | the Ensembl sources, the definitions, and the builder's validation records |
| `SHA256SUMS` | checksums of every file |

A component is a connected component of the Ensembl Compara 116 human-mouse orthology graph: the human and mouse
descendants of one ancestral gene. A one-to-one pair is a component with exactly one protein-coding gene in each
species. `provenance/ORTHOGROUP_DEFINITION.md` has the full definitions.

`mi/` reads the table from the path `ORTHOLOGS` in `mi/config.sh` (a byte-identical copy on the lab's storage) and
matches the one-to-one pairs to the objects' genes by Ensembl gene ID.

Sources: Ensembl release 116, retrieved on 7 August 2026; see `provenance/ENSEMBL_116_SOURCES.md`.
