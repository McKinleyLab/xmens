# orthology: the human-mouse ortholog table of record (Ensembl 116)

One copy of the human-mouse orthology for every script and notebook in the project, with a small Python module to
query it. The table is the one every cross-species analysis in the paper used. It is byte-identical to
`xmens_final/resources/orthology/` and `human_ivd_repro/orthogroups/`, and `check_orthology.py` confirms that.

| file | contents |
|---|---|
| `ensembl_human_mouse_components.tsv` | the table: one row per protein-coding human or mouse gene with a human-mouse homology call in Ensembl Compara 116 (36,334 genes in 17,030 components); sha256 `7c69c3a5…` |
| `ortholog_pairs_one2one.csv` | the 16,275 one-to-one symbol pairs (`human`, `mouse`) used in the paper; identical to the file rule `mi_overlap` writes in `xmens_final` |
| `ortholog_pairs_one2one_gene_ids.csv` | all 16,314 one-to-one pairs by Ensembl gene ID, with symbols |
| `orthology.py` | the module: `components`, `one_to_one`, `human_to_mouse`, `mouse_to_human`, `homology` |
| `check_orthology.py` | checks this folder against its checksums, the other copies and the IVD consensus mapping |
| `provenance/` | how the table was built: the Ensembl sources, the definitions, and the builder's validation records |
| `build/` | the builder and the frozen Ensembl 116 sources; `bash build/rebuild.sh` rebuilds the table byte for byte |
| `SHA256SUMS` | checksums of every file |

## Using it

From Python (a script, or a notebook's first code cell):

```python
import sys
sys.path.insert(0, "/n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/orthology")
import orthology

pairs = orthology.one_to_one()                          # DataFrame: human, mouse (16,275 rows)
mouse = orthology.human_to_mouse(adata.var_names)       # Series: the mouse symbol of each human gene, NaN if none
human = orthology.mouse_to_human(ids, by="gene_id")     # the same by Ensembl gene ID
orthology.homology(["PRL", "LEFTY2"])                   # every mouse homolog, for one-to-many genes too
```

The module needs only pandas; any of the project environments has it. It reads the table once per session and refuses a
table whose checksum differs from the table of record.

From a shell script or a `config.sh`, point at the files:

```bash
ORTHOLOGY=/n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/orthology
ORTHOLOG_PAIRS=$ORTHOLOGY/ortholog_pairs_one2one.csv
ORTHOLOG_COMPONENTS=$ORTHOLOGY/ensembl_human_mouse_components.tsv
```

From R: `read.csv(file.path(orthology, "ortholog_pairs_one2one.csv"))`, or `read.delim` for the table (use
`na.strings = ""`).

## Definitions

- **Component:** a connected component of the Ensembl 116 human-mouse orthology graph, meaning the human and mouse
  descendants of one ancestral gene. Of the 17,030 components, 16,314 are one-to-one, 202 have extra human copies, 397
  have extra mouse copies and 117 have both. `homology_type` is Ensembl's call for each gene, and `human_copy_number` and
  `mouse_copy_number` give the component's size. `provenance/ORTHOGROUP_DEFINITION.md` has the full definition.
- **One-to-one pair:** a component with exactly one gene in each species, all protein-coding.
  - By gene ID there are 16,314.
  - By symbol there are 16,275. 37 pairs drop out because a gene has no symbol. Then mouse `Tsc22d4` drops out on both
    sides, because it names two mouse genes partnered with TSC22D4 and with SPACDR.
  - This is the rule of `mi_overlap` (`xmens_final/workflow/scripts/mi/mi_overlap.py`), which defined the pairs in the
    paper.
- **Symbols naming two genes:** five mouse symbols name two Ensembl genes each: `4933427D14Rik`, `Ndor1`, `Nrg1`,
  `St6galnac2` and `Tsc22d4`. By symbol, the first four map to the gene whose human partner has a symbol. Match by
  gene ID where an object has them.
- **Not one-to-one:** genes such as PRL (1 human, 26 mouse genes) have no one-to-one partner. `homology()` lists all of
  their homologs.

## Who reads it

| user | reads | note |
|---|---|---|
| `xmens_final` rule `mi_overlap` (MI pseudotime overlap, Methods) | the table, through `orthologs:` in `config.yaml`, which points at `human_ivd_repro/orthogroups/` | writes `ortholog_pairs_one2one.csv` |
| `scenic` | a copy of `ortholog_pairs_one2one.csv` in `resources/`, and the table from `xmens_final/resources/orthology/` | both identical to the files here |
| IVD consensus mapping (`derive_consensus123_orthology.py`) | the table by Ensembl gene ID, plus OMA HOGs from `human_ivd_repro/orthogroups/` | agrees with `homology()` for all 123 genes |
| Stro-7 cosine job (`cosine_all_nuclei_job.py`, 29 Oct) | the table, building its own symbol lookup | gets 16,271 pairs, see below |
| other `xmens_final` code | the pairs file or the table | `regulon_rho_pseudotime.py`, the Methods-number scripts, `paper_code/ulrich_mapping/run_mapping_mouse*.py` |

The Stro-7 cosine job builds its lookup without first removing genes that have no symbol. Its uniqueness test therefore
also removes NRG1, NDOR1, ST6GALNAC2 and KIAA0753: each of their mouse symbols also names a second mouse gene, whose
human partner has no symbol. The job uses the lookup only to check that C7 has a one-to-one mouse ortholog, so no
reported value depends on it.

New code should call `orthology.one_to_one()` rather than derive pairs itself.

## Provenance and rebuilding

- **Sources:** Ensembl release 116, retrieved on 7 August 2026: the GTF gene rows for GRCh38 and GRCm39, the Compara
  `protein_default` homology dumps and the gene-tree content. The official FTP checksums were verified at download.
  `provenance/ENSEMBL_116_SOURCES.md` and `SOURCE_MANIFEST.tsv` record what was downloaded and from where.
- **Build:** `build/build_ensembl116_components.py` builds the table from those frozen files, and
  `build/validate_ensembl116_components.py` recounts it independently. Both come verbatim from
  `human_ivd_repro/orthogroup_work/scripts`; the builder gained `--source-dir` and `--out`.
- **Acquisition:** `build/freeze_ensembl116_sources.py` did the original downloads; it isn't needed to rebuild.
- **Rebuild:** `bash build/rebuild.sh` rebuilds the table in a scratch folder and checks it, and both validation
  records, byte for byte against the ones here.

## Checks

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/orthology
/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python check_orthology.py   # a few seconds
bash build/rebuild.sh                                                                  # under a minute
```

Nothing in this folder should change. If the project ever moves to another Ensembl release, build it into a new folder
next to this one.
