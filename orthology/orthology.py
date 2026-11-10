"""Human-mouse orthology of record for the X-Mens project: Ensembl Compara release 116.

The table has one row per protein-coding gene, human or mouse, with a human-mouse homology call in Ensembl 116. Genes
are grouped into components, the connected components of the orthology graph: a component holds the human and mouse
descendants of one ancestral gene, so it is one-to-one, one-to-many or many-to-many. README.md in this folder has the
definitions, the provenance and the checks.

From a script or notebook:

    import sys
    sys.path.insert(0, "/n/lab_storage/eddy_lab/users/nhilgert/xmens_scripts/orthology")
    import orthology

    orthology.one_to_one()                    # the 16,275 human-mouse symbol pairs of record
    orthology.human_to_mouse(["PRL", "C7"])   # one-to-one mouse partner of each gene, NaN where there is none
    orthology.homology(["PRL"])               # every mouse homolog, whatever the orthology type
"""
import hashlib
import io
from functools import lru_cache
from pathlib import Path

import pandas as pd

FOLDER = Path(__file__).resolve().parent
TABLE = FOLDER / "ensembl_human_mouse_components.tsv"
TABLE_SHA256 = "7c69c3a567f7e94dca9ef0a42be94da29ce196b9bb2bc8db0339339db654cdb7"
OTHER = {"human": "mouse", "mouse": "human"}


def components(path=TABLE, check=True):
    """The whole table, one row per gene.

    With check=True the file must be byte-identical to the table of record, so a stale or edited copy fails loudly
    instead of quietly changing results.
    """
    return _read(str(path), check).copy()


@lru_cache(maxsize=None)
def _read(path, check):
    data = Path(path).read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if check and digest != TABLE_SHA256:
        raise ValueError(f"{path} has sha256 {digest}, not that of the table of record ({TABLE_SHA256})")

    # Only empty fields count as missing, so no gene symbol can be read as a missing value.
    table = pd.read_csv(io.BytesIO(data), sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    return table.astype({"human_copy_number": int, "mouse_copy_number": int})


def one_to_one(by="symbol"):
    """Strict one-to-one pairs: the components with exactly one gene in each species.

    by="symbol": the 16,275 (human, mouse) symbol pairs of record, identical to ortholog_pairs_one2one.csv in this
        folder and to the file rule mi_overlap writes in xmens_final. Pairs where either gene has no symbol are left
        out first. A symbol still in two pairs after that (mouse Tsc22d4, partner of both TSC22D4 and SPACDR) is then
        left out on both sides, so every symbol has exactly one partner.
    by="gene_id": all 16,314 pairs of Ensembl gene IDs, with the genes' symbols where they have one.
    """
    if by not in ("symbol", "gene_id"):
        raise ValueError(f'by must be "symbol" or "gene_id", not {by!r}')
    return _one_to_one(by).copy()


@lru_cache(maxsize=None)
def _one_to_one(by):
    table = components()
    one2one = table[(table.homology_type == "ortholog_one2one") & (table.biotype == "protein_coding")]
    human = one2one[one2one.species == "human"].set_index("component_id")
    mouse = one2one[one2one.species == "mouse"].set_index("component_id")
    if human.index.has_duplicates or not human.index.equals(mouse.index):
        raise ValueError("a one-to-one component does not hold exactly one gene per species")

    pairs = pd.DataFrame({
        "human_gene_id": human.ensembl_gene_id,
        "mouse_gene_id": mouse.ensembl_gene_id,
        "human": human.gene_symbol,
        "mouse": mouse.gene_symbol,
    }).reset_index(drop=True)
    if by == "gene_id":
        return pairs

    named = pairs.dropna(subset=["human", "mouse"]).drop_duplicates(["human", "mouse"])
    unique = ~named.human.duplicated(keep=False) & ~named.mouse.duplicated(keep=False)
    return named.loc[unique, ["human", "mouse"]].reset_index(drop=True)


def human_to_mouse(genes, by="symbol"):
    """The one-to-one mouse partner of each human gene, NaN where there is none.

    genes are symbols (by="symbol") or Ensembl gene IDs (by="gene_id"); the result is indexed by them, in their order.
    """
    return _partner(genes, "human", by)


def mouse_to_human(genes, by="symbol"):
    """The one-to-one human partner of each mouse gene, NaN where there is none; see human_to_mouse."""
    return _partner(genes, "mouse", by)


def _partner(genes, species, by):
    pairs = one_to_one(by)
    suffix = "_gene_id" if by == "gene_id" else ""
    lookup = pairs.set_index(species + suffix)[OTHER[species] + suffix]
    genes = [genes] if isinstance(genes, str) else list(genes)
    return lookup.reindex(genes).rename_axis(species).rename(OTHER[species])


def homology(genes, species="human", by="symbol"):
    """Every homolog in the other species, whatever the orthology type (one-to-one, one-to-many, many-to-many).

    One row per query and matching gene. A query can match two genes (five mouse symbols, such as Nrg1, each name two
    Ensembl genes), and a query that is not in the table keeps one row with everything else empty. The homologs'
    symbols and gene IDs are joined with "; ", in table order.
    """
    if species not in OTHER:
        raise ValueError(f'species must be "human" or "mouse", not {species!r}')
    if by not in ("symbol", "gene_id"):
        raise ValueError(f'by must be "symbol" or "gene_id", not {by!r}')
    key = "gene_symbol" if by == "symbol" else "ensembl_gene_id"

    table = components()
    homologs = table[table.species == OTHER[species]].groupby("component_id", sort=False).agg(
        homolog_symbols=("gene_symbol", lambda symbols: "; ".join(symbols.dropna())),
        homolog_gene_ids=("ensembl_gene_id", "; ".join),
    )
    genes = [genes] if isinstance(genes, str) else list(genes)
    matches = table[(table.species == species) & table[key].isin(genes)].join(homologs, on="component_id")

    result = pd.DataFrame({"query": genes}).merge(matches, left_on="query", right_on=key, how="left")
    columns = ["query", "ensembl_gene_id", "gene_symbol", "component_id", "homology_type",
               "human_copy_number", "mouse_copy_number", "homolog_symbols", "homolog_gene_ids"]
    return result[columns].astype({"human_copy_number": "Int64", "mouse_copy_number": "Int64"})
