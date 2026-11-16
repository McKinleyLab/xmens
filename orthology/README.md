# orthology

The human-mouse ortholog table (Ensembl Compara 116).

| file | |
|---|---|
| `ensembl_human_mouse_components.tsv` | genes with a human-mouse homolog, grouped into components; read by `mi/` step 4 |
| `ortholog_pairs_one2one.csv` | the one-to-one pairs, by symbol |
| `orthology.py` | functions to query the table |
| `build/` | the code and frozen Ensembl 116 sources that make the table (`bash build/rebuild.sh`) |
| `provenance/` | the sources, the definitions and the builder's checks |
| `SHA256SUMS` | checksums |
