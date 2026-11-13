"""Step 5: the supplementary table of pseudotime MI for one-to-one human-mouse orthologs, SuppTable_pseudotime_MI.xlsx.

Sheets of the workbook:
  overlap          front sheet with live Excel formulas: change the percentile cutoff in B4, and the MI thresholds,
                   top-gene counts, N, K, n, k, expected overlap, fold enrichment and hypergeometric p recalculate from
                   the other sheets
  README           what each column means, and the input files with their sha256
  human_mouse_mi   one row per one-to-one pair with at least one gene tested; MI in bits
  human_genes, mouse_genes
                   every tested gene of the species, genome-wide; MI in bits
  p_value_terms    helper of the overlap sheet: the log hypergeometric probability of each overlap x, summed on the log
                   scale, as a p this small would underflow as 1 - CDF
N, K, n and k must agree with step 4's methods_overlap.json, or the script stops; the numbers in README come from the
inputs. Writes the folder --out-dir: the workbook, its human_mouse_mi, human_genes and mouse_genes sheets as CSVs, and a
provenance JSON (the inputs and their sha256).
Same tables as build_human_mouse_mi_table.py, which built the approved workbook of 3 Nov 2026 (record/5_table/), with
the code reorganized for reading. Since 12 Nov 2026 the workbook no longer carries the 2025 numbers beside the live
values (former column C of overlap, sheet checks); compare_with_record.py compares a run with the record.
"""
import argparse
import datetime
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.workbook.properties import CalcProperties
from scipy.stats import hypergeom

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--mi-mouse", required=True, help="step 3's mi_mouse.csv (MI in nats)")
parser.add_argument("--mi-human", required=True, help="step 3's mi_human.csv (MI in nats)")
parser.add_argument("--pairs", required=True, help="step 4's ortholog_pairs_one2one.csv")
parser.add_argument("--overlap", required=True, help="step 4's methods_overlap.json; N, K, n, k must agree with it")
parser.add_argument("--inputs", required=True, help="step 1's inputs folder (<species>/pseudotime.npy), for README")
parser.add_argument("--out-dir", required=True, help="output folder; must not exist yet")
args = parser.parse_args()

out_dir = args.out_dir.rstrip("/")
if os.path.exists(out_dir):
    sys.exit(f"{out_dir} exists; not overwriting it")

NAME = "SuppTable_pseudotime_MI"
SPECIES = ("human", "mouse")
LN2 = float(np.log(2))          # MI in bits = MI in nats / ln 2
Q_THRESH, PERCENTILE = 0.05, 90


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


# ---- inputs ----------------------------------------------------------------------------------------------------------
inputs = {"human": args.mi_human, "mouse": args.mi_mouse}
mi = {}
for species in SPECIES:
    mi[species] = pd.read_csv(inputs[species])
    if mi[species]["gene"].duplicated().any():
        sys.exit(f"duplicate genes in {inputs[species]}")
    mi[species] = mi[species].set_index("gene")
pairs = pd.read_csv(args.pairs)[["human", "mouse"]]
if pairs["human"].duplicated().any() or pairs["mouse"].duplicated().any():
    sys.exit(f"{args.pairs}: the pairs are not one-to-one")

# ---- per species: significant genes, the top-gene threshold, rank and percentile among significant genes --------------
stats = {}
for species in SPECIES:
    table = mi[species]
    significant = table[table["q"] < Q_THRESH]
    threshold = float(np.percentile(significant["mi"], PERCENTILE))          # nats
    rank = significant["mi"].rank(ascending=False, method="first")         # 1 = highest MI
    stats[species] = {
        "significant": significant,
        "threshold": threshold,
        "rank": rank,
        "percentile": 100 * (1 - (rank - 1) / len(significant)),
        "top": set(significant.index[significant["mi"] > threshold]),
    }


def label(in_human, in_mouse):
    """both / human / mouse / neither, per pair."""
    return np.where(in_human & in_mouse, "both",
                    np.where(in_human, "human", np.where(in_mouse, "mouse", "neither")))


# ---- the pair table (sheet human_mouse_mi) ---------------------------------------------------------------------------
pair_table = pairs.rename(columns={"human": "human_gene", "mouse": "mouse_gene"})
tested, significant, top = {}, {}, {}
for species in SPECIES:
    genes = pair_table[f"{species}_gene"]
    table, st = mi[species], stats[species]
    tested[species] = genes.isin(table.index).to_numpy()
    significant[species] = genes.isin(st["significant"].index).to_numpy()
    top[species] = genes.isin(st["top"]).to_numpy()
    pair_table[f"{species}_mi_rank"] = st["rank"].reindex(genes).to_numpy()
    pair_table[f"{species}_mi_percentile"] = st["percentile"].reindex(genes).to_numpy()
    pair_table[f"{species}_mi"] = table["mi"].reindex(genes).to_numpy() / LN2
    pair_table[f"{species}_q"] = table["q"].reindex(genes).to_numpy()
    pair_table[f"{species}_p"] = table["p"].reindex(genes).to_numpy()
pair_table["top_10_pct_mi"] = label(top["human"], top["mouse"])
pair_table["significant"] = label(significant["human"], significant["mouse"])
pair_table["expr_1_pct"] = label(tested["human"], tested["mouse"])
pair_table["ortholog_univ"] = significant["human"] & significant["mouse"]

# Pairs with neither gene tested are left out; rows sorted by top_10_pct_mi (both, human, mouse, neither), then by
# human rank, then mouse rank.
n_pairs = len(pair_table)
n_neither_tested = int((pair_table["expr_1_pct"] == "neither").sum())
pair_table = pair_table[pair_table["expr_1_pct"] != "neither"]
columns = ["human_gene", "mouse_gene", "top_10_pct_mi", "human_mi_rank", "mouse_mi_rank", "human_mi_percentile",
           "mouse_mi_percentile", "human_mi", "mouse_mi", "human_q", "mouse_q", "human_p", "mouse_p", "significant",
           "expr_1_pct", "ortholog_univ"]
TOP_ORDER = {"both": 0, "human": 1, "mouse": 2, "neither": 3}
pair_table = (pair_table.assign(top_order=pair_table["top_10_pct_mi"].map(TOP_ORDER))
              .sort_values(["top_order", "human_mi_rank", "mouse_mi_rank", "human_gene"], na_position="last"))
pair_table = pair_table[columns].reset_index(drop=True)
for species in SPECIES:
    pair_table[f"{species}_mi_rank"] = pair_table[f"{species}_mi_rank"].astype("Int64")

# ---- the gene tables (sheets human_genes and mouse_genes): every tested gene, by MI --------------------------------------
gene_tables = {}
for species in SPECIES:
    genes = mi[species].reset_index().sort_values("mi", ascending=False, kind="stable").reset_index(drop=True)
    genes["significant"] = genes["q"] < Q_THRESH
    genes["rank_among_significant"] = stats[species]["rank"].reindex(genes["gene"]).astype("Int64").to_numpy()
    genes["mi"] = genes["mi"] / LN2
    genes["mi_if_significant"] = genes["mi"].where(genes["significant"])
    gene_tables[species] = genes[["gene", "mi", "p", "q", "significant", "rank_among_significant", "mi_if_significant"]]

# ---- the overlap test, which must agree with step 4 ------------------------------------------------------------------
background = pair_table[pair_table["ortholog_univ"]]
N = len(background)
K = int(background["top_10_pct_mi"].isin(["mouse", "both"]).sum())
n = int(background["top_10_pct_mi"].isin(["human", "both"]).sum())
k = int((background["top_10_pct_mi"] == "both").sum())
log10_p = float(hypergeom.logsf(k - 1, N, K, n) / np.log(10))
mouse, human = stats["mouse"], stats["human"]
with open(args.overlap) as f:
    step4 = json.load(f)
if (step4["N"], step4["K"], step4["n"], step4["k"]) != (N, K, n, k):
    sys.exit(f"N, K, n, k = {N}, {K}, {n}, {k} here, but {step4['N']}, {step4['K']}, {step4['n']}, {step4['k']} in "
             f"{args.overlap}")

# ---- nuclei with a pseudotime: the rows of step 1's matrices ---------------------------------------------------------
n_nuclei = {species: len(np.load(os.path.join(args.inputs, species, "pseudotime.npy"))) for species in SPECIES}

# ---- README ---------------------------------------------------------------------------------------------------------
built = datetime.date.today().isoformat()
tested_in = pair_table["expr_1_pct"].value_counts()
readme = pd.DataFrame([
    ("Table", "Pseudotime-associated genes (mutual information with pseudotime) in human and mouse fibroblasts and "
              "decidual cells, for strict one-to-one human-mouse orthologs. New version of human_mouse_mi (2 Jul 2026), "
              "rebuilt on the final pseudotime run."),
    ("Built", f"{built} by xmens_mi step 5 (scripts/supp_table.py) from the inputs below (with their sha256)."),
    ("Pseudotime", "obs['pseudotime'] of the final fibroblast/decidual objects (rank-normalized diffusion pseudotime): "
                   f"{n_nuclei['mouse']:,} mouse and {n_nuclei['human']:,} human nuclei. Apoptotic decidual nuclei "
                   "have no pseudotime and are not used."),
    ("MI", "Per gene, MI between log-normalized expression and pseudotime (sklearn mutual_info_regression, n_neighbors = "
           "3, fixed random state); p from 1,000 pseudotime permutations, (count >= observed + 1) / 1,001; BH-adjusted "
           "q. MI is in bits."),
    ("Rows", f"Ensembl Compara 116 one-to-one pairs ({n_pairs:,}) with at least one gene tested: {len(pair_table):,} "
             f"pairs (tested in both {tested_in.get('both', 0):,}, human only {tested_in.get('human', 0):,}, mouse only "
             f"{tested_in.get('mouse', 0):,}); {n_neither_tested:,} pairs with neither gene tested are omitted. Sorted by "
             "top_10_pct_mi (both, human, mouse, neither), then human rank, then mouse rank."),
    ("Untested genes", "A gene detected in < 1% of nuclei was not tested and has no MI, p, q, rank or percentile."),
    ("human_gene / mouse_gene", "Ensembl Compara 116 one-to-one ortholog pair (homology type ortholog_one2one, "
                                "protein-coding genes), matched to the objects' genes by Ensembl gene ID; gene "
                                "names as in the objects."),
    ("top_10_pct_mi", "Species in which the gene is a top pseudotime-associated gene: MI above the 90th percentile of MI "
                      f"among that species' significant genes (mouse {mouse['threshold'] / LN2:.3f} bits; human "
                      f"{human['threshold'] / LN2:.3f} bits). both / human / mouse / neither. Genome-wide sets "
                      f"({len(mouse['top']):,} mouse, {len(human['top']):,} human); within the background they give K, "
                      "n and the overlap k."),
    ("human_mi_rank / mouse_mi_rank", "Rank of the gene's MI among all significant genes of its species, genome-wide "
                                      "(1 = highest). Empty if not significant."),
    ("human_mi_percentile / mouse_mi_percentile", "100 x (1 - (rank - 1) / number of significant genes in the "
                                                  "species). A top gene has percentile >= 90; the gene exactly at the "
                                                  "90th-percentile MI is not counted as top (strict >)."),
    ("human_mi / mouse_mi", "Mutual information with pseudotime (bits)."),
    ("human_q / mouse_q", "Benjamini-Hochberg q over all tested genes of the species."),
    ("human_p / mouse_p", "Permutation p-value (1,000 permutations; minimum 1/1,001)."),
    ("significant", "Species in which the gene is tested and q < 0.05: both / human / mouse / neither."),
    ("expr_1_pct", "Species in which the gene was detected in >= 1% of nuclei and so tested: both / human / mouse."),
    ("ortholog_univ", "TRUE if both genes are significant: the background of the cross-species comparison "
                      f"(N = {N:,} pairs; hypergeometric test on the overlap sheet)."),
    ("Changes from the 2 Jul version", "Final pseudotime and MI of record; ortholog pairs from Ensembl Compara 116 (the "
                                       "July table had 15,134 pairs from the earlier homology set, including "
                                       "non-coding genes such as NEAT1 and MALAT1); untested genes now "
                                       "have no MI (the July run scored every gene and flagged expression afterwards); "
                                       "'significant' is now q < 0.05 among tested genes (the July file repeated "
                                       "expr_1_pct in this column); top_10_pct_mi is MI above the 90th percentile "
                                       "among significant genes."),
    ("Input: MI, mouse", f"{args.mi_mouse} (sha256 {sha256(args.mi_mouse)})"),
    ("Input: MI, human", f"{args.mi_human} (sha256 {sha256(args.mi_human)})"),
    ("Input: ortholog pairs", f"{args.pairs} (sha256 {sha256(args.pairs)})"),
    ("Sheet: overlap", "Front sheet. Live Excel formulas: change the cutoff in B4 and the thresholds, top-gene counts, "
                       "background N, K, n, overlap k, expected overlap, fold enrichment and the hypergeometric p "
                       "recalculate from the other sheets."),
    ("Sheet: human_genes / mouse_genes", "Every tested gene of the species (genome-wide, not only orthologs): MI (bits), "
                                         "p, q, significant (q < 0.05), rank among significant genes, and "
                                         "mi_if_significant (the column the percentile threshold is taken over)."),
    ("Sheet: p_value_terms", "Helper for the overlap sheet: log hypergeometric probability of each overlap x (log-gamma "
                             "form), so the upper tail P(X >= k) is summed without underflow (a p this small cannot "
                             "be computed as 1 - CDF)."),
], columns=["field", "description"])

# ---- the workbook -----------------------------------------------------------------------------------------------------
partial = os.path.join(os.path.dirname(out_dir), ".partial." + os.path.basename(out_dir))
os.makedirs(partial, exist_ok=True)
bold = Font(bold=True)
with pd.ExcelWriter(os.path.join(partial, f"{NAME}.xlsx"), engine="openpyxl") as writer:
    readme.to_excel(writer, sheet_name="README", index=False)
    pair_table.to_excel(writer, sheet_name="human_mouse_mi", index=False)
    gene_tables["human"].to_excel(writer, sheet_name="human_genes", index=False)
    gene_tables["mouse"].to_excel(writer, sheet_name="mouse_genes", index=False)
    wb = writer.book

    # Column widths, bold headers, frozen header row, and number formats of the pair table.
    widths = {
        "README": {"A": 34, "B": 140},
        "human_mouse_mi": {column: 15 for column in "ABCDEFGHIJKLMNOP"},
        "human_genes": {column: 16 for column in "ABCDEFG"},
        "mouse_genes": {column: 16 for column in "ABCDEFG"},
    }
    for sheet, sheet_widths in widths.items():
        ws = wb[sheet]
        for column, width in sheet_widths.items():
            ws.column_dimensions[column].width = width
        for cell in ws[1]:
            cell.font = bold
        ws.freeze_panes = "A2"
    for row in wb["README"].iter_rows(min_row=2):
        row[0].alignment = Alignment(vertical="top")
        row[1].alignment = Alignment(wrap_text=True, vertical="top")
    formats = {"F": "0.00", "G": "0.00", "H": "0.0000", "I": "0.0000",
               "J": "0.00E+00", "K": "0.00E+00", "L": "0.00E+00", "M": "0.00E+00"}
    for column, number_format in formats.items():
        for cell in wb["human_mouse_mi"][column][1:]:
            cell.number_format = number_format

    # Sheet p_value_terms: log P(X = x) of the hypergeometric overlap for x = 0 to N + 50, using the overlap sheet's
    # N, K, n and k (log binomial coefficients from GAMMALN); F2 sums the upper tail x >= k as a log-sum-exp.
    terms = wb.create_sheet("p_value_terms")
    cell_N, cell_K, cell_n, cell_k = "overlap!$B$17", "overlap!$B$18", "overlap!$B$19", "overlap!$B$20"
    x_max = int(pair_table["ortholog_univ"].sum()) + 50

    def log_choose(a, b):
        return f"(GAMMALN({a}+1)-GAMMALN({b}+1)-GAMMALN({a}-{b}+1))"

    terms.append(["x", "log P(X = x), x in the upper tail [k, min(K, n)]", "exp(log P - max)", None, "max log P",
                  "log10 p = log10 P(X >= k)"])
    for x_value in range(x_max + 1):
        row = x_value + 2
        x = f"A{row}"
        log_pmf = (f"{log_choose(cell_K, x)}+{log_choose(f'({cell_N}-{cell_K})', f'({cell_n}-{x})')}"
                   f"-{log_choose(cell_N, cell_n)}")
        in_tail = f"AND({x}>={cell_k},{x}<=MIN({cell_K},{cell_n}),{cell_n}-{x}<={cell_N}-{cell_K})"
        terms.append([x_value, f"=IF({in_tail},{log_pmf},-1E+300)", f"=EXP(B{row}-$E$2)"])
    terms["E2"] = f"=MAX($B$2:$B${x_max + 2})"
    terms["F2"] = f"=($E$2+LN(SUM($C$2:$C${x_max + 2})))/LN(10)"
    for cell in terms[1]:
        cell.font = bold
    for column, width in {"B": 44, "C": 18, "E": 14, "F": 26}.items():
        terms.column_dimensions[column].width = width

    # Sheet overlap, the live front sheet: each row is (row, quantity, formula, number format).
    overlap = wb.create_sheet("overlap", 0)
    last_pair, last_human, last_mouse = len(pair_table) + 1, len(gene_tables["human"]) + 1, len(gene_tables["mouse"]) + 1
    univ = f"human_mouse_mi!$P$2:$P${last_pair}"            # ortholog_univ
    mi_h = f"human_mouse_mi!$H$2:$H${last_pair}"            # human_mi
    mi_m = f"human_mouse_mi!$I$2:$I${last_pair}"            # mouse_mi
    sig_h = f"human_genes!$G$2:$G${last_human}"             # mi_if_significant
    sig_m = f"mouse_genes!$G$2:$G${last_mouse}"
    genes_h = f"human_genes!$A$2:$A${last_human}"
    genes_m = f"mouse_genes!$A$2:$A${last_mouse}"
    overlap["A1"] = ("Cross-species overlap of top pseudotime-associated genes (live: computed from the other sheets of "
                     "this workbook)")
    overlap["A2"] = "Change the cutoff in B4; every value below recalculates from the other sheets."
    rows = [
        (4, "Top-gene cutoff: percentile of MI among each species' significant genes", 90, "0"),
        (6, "Mouse MI threshold (bits)", f"=PERCENTILE({sig_m},$B$4/100)", "0.000"),
        (7, "Human MI threshold (bits)", f"=PERCENTILE({sig_h},$B$4/100)", "0.000"),
        (10, "Mouse genes tested", f"=COUNTA({genes_m})", "#,##0"),
        (11, "Human genes tested", f"=COUNTA({genes_h})", "#,##0"),
        (12, "Mouse significant genes (q < 0.05)", f"=COUNT({sig_m})", "#,##0"),
        (13, "Human significant genes (q < 0.05)", f"=COUNT({sig_h})", "#,##0"),
        (14, "Mouse top genes, genome-wide (MI > threshold)", f"=SUMPRODUCT(({sig_m}>B6)*1)", "#,##0"),
        (15, "Human top genes, genome-wide (MI > threshold)", f"=SUMPRODUCT(({sig_h}>B7)*1)", "#,##0"),
        (17, "Background pairs N (both genes significant)", f"=SUMPRODUCT(({univ}=TRUE)*1)", "#,##0"),
        (18, "Mouse top genes in the background, K", f"=SUMPRODUCT(({univ}=TRUE)*({mi_m}>B6))", "#,##0"),
        (19, "Human top genes in the background, n", f"=SUMPRODUCT(({univ}=TRUE)*({mi_h}>B7))", "#,##0"),
        (20, "Overlap k (top in both species)", f"=SUMPRODUCT(({univ}=TRUE)*({mi_m}>B6)*({mi_h}>B7))", "#,##0"),
        (21, "Expected overlap, K x n / N", "=B18*B19/B17", "0.0"),
        (22, "Expected overlap, % of the human top set", "=B21/B19", "0%"),
        (23, "Observed overlap, % of the human top set", "=B20/B19", "0%"),
        (24, "Fold enrichment, k / expected", "=B20/B21", "0.00"),
        (25, "log10 p, hypergeometric P(X >= k)", "=p_value_terms!$F$2", "0.0000"),
        (26, "Hypergeometric p, P(X >= k): full-precision number (reads 0 below about 1E-308)", "=10^B25",
         "0.00E+00"),
        (27, "Hypergeometric p, as text", '=TEXT(10^(B25-INT(B25)),"0.00")&" x 10^"&INT(B25)', "@"),
    ]
    overlap["A3"], overlap["B3"] = "quantity", "value"
    for row, quantity, formula, number_format in rows:
        overlap[f"A{row}"] = quantity
        overlap[f"B{row}"] = formula
        overlap[f"B{row}"].number_format = number_format
    overlap["A29"] = ("N, K, n and k are counted in human_mouse_mi; the thresholds and genome-wide counts in human_genes / "
                      "mouse_genes; p is summed in p_value_terms. MI thresholds use PERCENTILE (linear interpolation, as "
                      "numpy.percentile).")
    for cell in overlap[3]:
        cell.font = bold
    overlap["A1"].font = Font(bold=True, size=12)
    overlap["B4"].fill = PatternFill("solid", fgColor="FFF2CC")
    overlap["B4"].font = bold
    for row in range(17, 28):            # the results of the overlap test, in bold: N, K, n, k, expected, %, fold, p
        overlap[f"B{row}"].font = bold
    for column, width in {"A": 84, "B": 18}.items():
        overlap.column_dimensions[column].width = width

    # Sheet order, the front sheet first; Excel recalculates every formula when the file is opened.
    wb._sheets = [wb[sheet] for sheet in
                  ("overlap", "README", "human_mouse_mi", "human_genes", "mouse_genes", "p_value_terms")]
    wb.active = 0
    wb.calculation = CalcProperties(fullCalcOnLoad=True)

# ---- the CSVs and the provenance record; then the folder is renamed when complete -------------------------------------
pair_table.to_csv(os.path.join(partial, f"{NAME}__human_mouse_mi.csv"), index=False)
for species in SPECIES:
    gene_tables[species].to_csv(os.path.join(partial, f"{NAME}__{species}_genes.csv"), index=False)
provenance = {
    "built": built,
    "inputs": {key: {"path": path, "sha256": sha256(path)}
               for key, path in (("mi_mouse", args.mi_mouse), ("mi_human", args.mi_human), ("pairs", args.pairs))},
    "rows": len(pair_table),
    "omitted_neither_tested": n_neither_tested,
    "N": N, "K": K, "n": n, "k": k,
    "log10_p": log10_p,
    "N_K_n_k_agree_with_step4": True,
    "nuclei_with_pseudotime": n_nuclei,
}
with open(os.path.join(partial, f"{NAME}_provenance.json"), "w") as f:
    json.dump(provenance, f, indent=1)
os.rename(partial, out_dir)
print(f"N {N:,}, K {K:,}, n {n:,}, k {k}, as in step 4; log10 p {log10_p:.4f}; nuclei with a pseudotime: mouse "
      f"{n_nuclei['mouse']:,}, human {n_nuclei['human']:,}")
print(f"{len(pair_table):,} pairs ({n_neither_tested:,} with neither gene tested left out); wrote {out_dir}", flush=True)
