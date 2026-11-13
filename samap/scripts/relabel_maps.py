"""Step 0c: the SAMap homology maps, with each gene under the name it has in the objects SAMap reads.

Step 0b names the genes of the maps by the gene_symbol of the Ensembl release 113 proteomes. The objects name genes as
the annotation they were counted with does, which is older: a gene renamed since (mouse Sept2, now Septin2; human
C1orf112, now FIRRM) is in the maps under its new name only, so SAMap finds no homology for it. This step renames such
genes in both maps to the objects' names, matched by Ensembl gene ID. Only the first two fields of a line (the two gene
names) can change; the rest of every line is kept byte for byte.

A map name S becomes the object's name N when
  - of the Ensembl 113 genes named S (gene IDs from the proteome FASTA headers), exactly one is in the object, under N;
  - N differs from S, no other map name becomes N, and N is not a map name that stays (two genes would merge).
Writes --out, which must not exist: hsmm/hs_to_mm.txt and hsmm/mm_to_hs.txt (the relabeled maps), renamed_genes.tsv
(one row per renamed gene), not_renamed.tsv (candidates left as they were, with the reason) and relabel.json (the
inputs with their sha256, and the counts).
"""
import argparse
import datetime
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict

import h5py

SPECIES = ("human", "mouse")
# The map files, and which species each of a line's first two fields belongs to.
MAP_FILES = {"hs_to_mm.txt": ("human", "mouse"), "mm_to_hs.txt": ("mouse", "human")}

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--maps", required=True,
                    help="folder holding hsmm/hs_to_mm.txt and hsmm/mm_to_hs.txt (step 0b)")
parser.add_argument("--human-proteome", required=True, help="the human proteome FASTA the maps were made from")
parser.add_argument("--mouse-proteome", required=True, help="the mouse proteome FASTA the maps were made from")
parser.add_argument("--human-proteome-sha256", default="", help="if given, the human proteome must have this sha256")
parser.add_argument("--mouse-proteome-sha256", default="", help="if given, the mouse proteome must have this sha256")
parser.add_argument("--human-h5ad", required=True,
                    help="the human object SAMap reads (step 1's hs_fib_dec_all_raw.h5ad)")
parser.add_argument("--mouse-h5ad", required=True,
                    help="the mouse object SAMap reads (step 1's mm_fib_dec_all_raw.h5ad)")
parser.add_argument("--gene-id-key", default="gene_ids",
                    help="var column with the Ensembl gene IDs (default: gene_ids)")
parser.add_argument("--out", required=True, help="output folder; must not exist yet")
args = parser.parse_args()

out = args.out.rstrip("/")
if os.path.exists(out):
    sys.exit(f"{out} exists; not overwriting it")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def proteome_genes(path):
    """Ensembl gene ID (unversioned) -> gene_symbol, from the FASTA headers; proteins without a symbol are skipped."""
    symbols = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                gene = re.search(r" gene:(\S+)", line)
                symbol = re.search(r" gene_symbol:(\S+)", line)
                if gene and symbol:
                    symbols[gene.group(1).split(".")[0]] = symbol.group(1)
    return symbols


def object_genes(path):
    """Ensembl gene ID (unversioned) -> the object's gene name, for the IDs that occur once in the object."""

    def strings(item):
        if isinstance(item, h5py.Group):  # a categorical column
            return item["categories"][:].astype(str)[item["codes"][:]]
        return item[:].astype(str)

    with h5py.File(path, "r") as handle:
        var = handle["var"]
        names = strings(var[var.attrs.get("_index", "_index")])
        ids = [value.split(".")[0] for value in strings(var[args.gene_id_key])]
    counts = Counter(ids)
    return {gene_id: name for gene_id, name in zip(ids, names) if counts[gene_id] == 1}


def map_names(paths):
    """The gene names of each species in the two maps."""
    names = {species: set() for species in SPECIES}
    for file_name, (first, second) in MAP_FILES.items():
        with open(paths[file_name]) as handle:
            for line in handle:
                a, b, _ = line.split("\t", 2)
                names[first].add(a)
                names[second].add(b)
    return names


def renames(names, proteome, genes):
    """{map name: object name} for one species, and {map name: (object name, reason)} for the candidates left out."""
    ids_by_symbol = defaultdict(list)
    for gene_id, symbol in proteome.items():
        ids_by_symbol[symbol].append(gene_id)
    proposed, left = {}, {}
    for name in names:
        in_object = {genes[gene_id] for gene_id in ids_by_symbol.get(name, []) if gene_id in genes}
        if len(in_object) > 1:
            left[name] = (";".join(sorted(in_object)), "several genes of this name in the object")
        elif in_object and (target := in_object.pop()) != name:
            proposed[name] = target
    # Drop a renaming whose target another map name also gets, or that is a map name which stays; repeat, as each drop
    # leaves one more name in place.
    while True:
        targets = Counter(proposed.values())
        staying = names - set(proposed)
        clash = {name: target for name, target in proposed.items() if targets[target] > 1 or target in staying}
        if not clash:
            break
        for name, target in clash.items():
            reason = "two map names would get this name" if targets[target] > 1 else "the name is taken in the map"
            left[name] = (target, reason)
            del proposed[name]
    final = [proposed.get(name, name) for name in names]
    if len(final) != len(set(final)):
        raise RuntimeError("relabeling would merge genes")
    return proposed, left


paths = {file_name: os.path.join(args.maps, "hsmm", file_name) for file_name in MAP_FILES}
proteomes = {"human": args.human_proteome, "mouse": args.mouse_proteome}
pins = {"human": args.human_proteome_sha256, "mouse": args.mouse_proteome_sha256}
h5ads = {"human": args.human_h5ad, "mouse": args.mouse_h5ad}
for species in SPECIES:
    if pins[species] and sha256(proteomes[species]) != pins[species]:
        sys.exit(f"{proteomes[species]}: sha256 differs from the pinned {pins[species]}")

names = map_names(paths)
renamed, left, counts = {}, {}, {}
for species in SPECIES:
    renamed[species], left[species] = renames(names[species], proteome_genes(proteomes[species]),
                                              object_genes(h5ads[species]))
    counts[species] = {"map_names": len(names[species]), "renamed": len(renamed[species]),
                       "left_as_they_were": len(left[species])}
    print(f"{species}: {len(names[species]):,} map names; {len(renamed[species]):,} renamed to the object's name; "
          f"{len(left[species]):,} candidates left as they were", flush=True)

# Written into a temporary folder that is renamed when complete, so an interrupted run never leaves a folder that
# looks finished.
partial = os.path.join(os.path.dirname(out), ".partial." + os.path.basename(out))
os.makedirs(os.path.join(partial, "hsmm"), exist_ok=True)
for file_name, (first, second) in MAP_FILES.items():
    changed = 0
    with open(paths[file_name]) as source, open(os.path.join(partial, "hsmm", file_name), "w") as target:
        for line in source:
            a, b, rest = line.split("\t", 2)
            new_a, new_b = renamed[first].get(a, a), renamed[second].get(b, b)
            changed += (new_a, new_b) != (a, b)
            target.write(f"{new_a}\t{new_b}\t{rest}")
    counts[file_name] = {"lines_changed": changed}
    print(f"{file_name}: {changed:,} lines with a renamed gene", flush=True)

with open(os.path.join(partial, "renamed_genes.tsv"), "w") as handle:
    handle.write("species\tmap_name_ensembl113\tobject_name\n")
    for species in SPECIES:
        for name, target in sorted(renamed[species].items()):
            handle.write(f"{species}\t{name}\t{target}\n")
with open(os.path.join(partial, "not_renamed.tsv"), "w") as handle:
    handle.write("species\tmap_name_ensembl113\tobject_name\treason\n")
    for species in SPECIES:
        for name, (target, reason) in sorted(left[species].items()):
            handle.write(f"{species}\t{name}\t{target}\t{reason}\n")
record = {
    "created": datetime.datetime.now().isoformat(timespec="seconds"),
    "script": {"path": os.path.abspath(__file__), "sha256": sha256(__file__)},
    "inputs": {key: {"path": path, "sha256": sha256(path)} for key, path in
               [*((f"maps/{name}", path) for name, path in paths.items()),
                *((f"proteome/{species}", path) for species, path in proteomes.items()),
                *((f"h5ad/{species}", path) for species, path in h5ads.items())]},
    "counts": counts,
}
with open(os.path.join(partial, "relabel.json"), "w") as handle:
    handle.write(json.dumps(record, indent=1) + "\n")
os.rename(partial, out)
print(f"wrote {out}", flush=True)
