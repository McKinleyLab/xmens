"""Checks this folder: checksums, the pair files, and agreement with the other copies and uses of the table.

Run after installing or changing anything here (a few seconds):
    /n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python check_orthology.py
Prints one line per check and exits with status 1 if any check fails. Files elsewhere on the cluster that are not
found are reported as skipped, so the script also runs off the cluster.
"""
import hashlib
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import orthology  # noqa: E402  (the module in this folder, whatever the working directory)

NH = Path("/n/lab_storage/eddy_lab/users/nhilgert")
TABLE_COPIES = [
    NH / "xmens_final/resources/orthology/ensembl_human_mouse_components.tsv",
    NH / "human_ivd_repro/orthogroups/ensembl_human_mouse_components.tsv",
]
PAIR_COPIES = [
    NH / "xmens_final/results/mi_dpt/overlap_of_record/ortholog_pairs_one2one.csv",
    NH / "xmens_final/results/mi_dpt/overlap/ortholog_pairs_one2one.csv",
    NH / "xmens_scripts/scenic/resources/ortholog_pairs_one2one.csv",
]
# The IVD consensus derivation maps its 123 genes by Ensembl gene ID; it is an independent use of the same table.
CONSENSUS = NH / "xmens_final/resources/ivd/expected/consensus/consensus123_stadtmauer_authors_criteria_mouse_orthology.csv"

# Counts from provenance/ensembl_116_component_validation.json, and the symbol pairs rule mi_overlap reports.
EXPECTED = dict(rows=36334, human_genes=17688, mouse_genes=18646, components=17030, one_to_one=16314,
                human_expansion=202, mouse_expansion=397, both_expanded=117, symbol_pairs=16275)

failures = []


def report(ok, name, detail=""):
    status = {True: "PASS", False: "FAIL", None: "SKIP"}[ok]
    print(f"{status}  {name}" + (f": {detail}" if detail else ""))
    if ok is False:
        failures.append(name)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


listed = [line.split(maxsplit=1) for line in (HERE / "SHA256SUMS").read_text().splitlines() if line.strip()]
changed = [name for digest, name in listed if sha256(HERE / name) != digest]
report(not changed, f"{len(listed)} files match SHA256SUMS", ", ".join(changed))

table = orthology.components()  # raises if the table is not byte-identical to the table of record
per_component = table.groupby("component_id")[["human_copy_number", "mouse_copy_number"]].first()
human_n, mouse_n = per_component.human_copy_number, per_component.mouse_copy_number
counts = dict(
    rows=len(table),
    human_genes=int((table.species == "human").sum()),
    mouse_genes=int((table.species == "mouse").sum()),
    components=len(per_component),
    one_to_one=int(((human_n == 1) & (mouse_n == 1)).sum()),
    human_expansion=int(((human_n > 1) & (mouse_n == 1)).sum()),
    mouse_expansion=int(((human_n == 1) & (mouse_n > 1)).sum()),
    both_expanded=int(((human_n > 1) & (mouse_n > 1)).sum()),
    symbol_pairs=len(orthology.one_to_one()),
)
report(counts == EXPECTED, "table counts as recorded", "" if counts == EXPECTED else str(counts))

for by, name in (("symbol", "ortholog_pairs_one2one.csv"), ("gene_id", "ortholog_pairs_one2one_gene_ids.csv")):
    regenerated = orthology.one_to_one(by).to_csv(index=False).encode()
    report(regenerated == (HERE / name).read_bytes(), f"{name} regenerates byte for byte")

for copies, here in ((TABLE_COPIES, HERE / "ensembl_human_mouse_components.tsv"),
                     (PAIR_COPIES, HERE / "ortholog_pairs_one2one.csv")):
    for copy in copies:
        if copy.exists():
            report(sha256(copy) == sha256(here), f"same bytes as {copy}")
        else:
            report(None, f"{copy}", "not found")

if CONSENSUS.exists():
    consensus = pd.read_csv(CONSENSUS, keep_default_na=False, na_values=[""])
    found = orthology.homology(consensus.human_ensembl_gene_id, by="gene_id").set_index("query")
    found = found.loc[consensus.human_ensembl_gene_id]
    same_type = found.homology_type.fillna("").to_numpy() == consensus.ensembl116_homology_type.fillna("").to_numpy()

    # The consensus table joins mouse members with ";" in alphabetical order; homology() joins with "; " in table order.
    def members(text, separator):
        return frozenset(text.split(separator)) if isinstance(text, str) and text else frozenset()
    same_members = [members(a, ";") == members(b, "; ") for a, b in
                    zip(consensus.ensembl116_mouse_members, found.homolog_symbols)]

    partner = orthology.human_to_mouse(consensus.human_ensembl_gene_id, by="gene_id").to_numpy()
    strict = consensus.strict_one_to_one.astype(bool).to_numpy()
    same_partner = (pd.notna(partner) == strict) & (
        ~strict | (partner == consensus.mouse_ensembl_gene_id.to_numpy()))
    detail = f"type {int(same_type.sum())}, members {sum(same_members)}, one-to-one partner {int(same_partner.sum())} of {len(consensus)}"
    report(bool(same_type.all() and all(same_members) and same_partner.all()),
           "agrees with the IVD consensus orthology table", detail)
else:
    report(None, f"{CONSENSUS}", "not found")

print(f"{len(failures)} failed")
sys.exit(1 if failures else 0)
