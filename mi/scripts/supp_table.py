"""Step 5: the supplementary table of pseudotime MI for one-to-one human-mouse orthologs, SuppTable_pseudotime_MI.xlsx.

Sheets of the workbook:
  overlap          front sheet with live Excel formulas: change the percentile cutoff in B4, and the MI thresholds,
                   top-gene counts, N, K, n, k, expected overlap, fold enrichment and hypergeometric p recalculate from
                   the other sheets; at 90 they reproduce the Methods (column C)
  README           what each column means, and the input files with their sha256
  human_mouse_mi   one row per one-to-one pair with at least one gene tested; MI in bits
  human_genes, mouse_genes
                   every tested gene of the species, genome-wide; MI in bits
  p_value_terms    helper of the overlap sheet: the log hypergeometric probability of each overlap x, summed on the log
                   scale, as p near 1e-84 would underflow as 1 - CDF
  checks           the Methods numbers recomputed here, and N, K, n, k against step 4's methods_overlap.json
Writes the folder --out-dir: the workbook, its human_mouse_mi, human_genes, mouse_genes and checks sheets as CSVs, and a
provenance JSON (the inputs and their sha256).
Same tables as build_human_mouse_mi_table.py, which built the approved workbook of 3 Nov 2026 (record/5_table/), with
the code reorganized for reading.
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
parser.add_argument("--overlap", required=True, help="step 4's methods_overlap.json, for the checks sheet")
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

# Pairs with neither gene tested are left out; rows sorted by human rank, then mouse rank.
n_pairs = len(pair_table)
n_neither_tested = int((pair_table["expr_1_pct"] == "neither").sum())
pair_table = pair_table[pair_table["expr_1_pct"] != "neither"]
columns = ["human_gene", "mouse_gene", "top_10_pct_mi", "human_mi_rank", "mouse_mi_rank", "human_mi_percentile",
           "mouse_mi_percentile", "human_mi", "mouse_mi", "human_q", "mouse_q", "human_p", "mouse_p", "significant",
           "expr_1_pct", "ortholog_univ"]
pair_table = pair_table.sort_values(["human_mi_rank", "mouse_mi_rank", "human_gene"], na_position="last")
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

# ---- checks against the numbers in the Methods --------------------------------------------------------------------------
background = pair_table[pair_table["ortholog_univ"]]
N = len(background)
K = int(background["top_10_pct_mi"].isin(["mouse", "both"]).sum())
n = int(background["top_10_pct_mi"].isin(["human", "both"]).sum())
k = int((background["top_10_pct_mi"] == "both").sum())
log10_p = float(hypergeom.logsf(k - 1, N, K, n) / np.log(10))
exponent = int(np.floor(log10_p))
mouse, human = stats["mouse"], stats["human"]
checks = pd.DataFrame([
    ("one-to-one pairs (Ensembl Compara 116)", "16,275", f"{n_pairs:,}"),
    ("genes tested, mouse / human", "14,261 / 15,008", f"{len(mi['mouse']):,} / {len(mi['human']):,}"),
    ("significant (q < 0.05), mouse / human", "14,261 / 14,887",
     f"{len(mouse['significant']):,} / {len(human['significant']):,}"),
    ("top genes genome-wide, mouse / human", "1,426 / 1,489", f"{len(mouse['top']):,} / {len(human['top']):,}"),
    ("top-gene MI threshold (bits), mouse / human", "0.102 / 0.039",
     f"{mouse['threshold'] / np.log(2):.3f} / {human['threshold'] / np.log(2):.3f}"),
    ("background pairs N (ortholog_univ = TRUE)", "9,587", f"{N:,}"),
    ("top mouse genes in background, K", "1,133", f"{K:,}"),
    ("top human genes in background, n", "1,201", f"{n:,}"),
    ("overlap k (top_10_pct_mi = both, in background)", "373", f"{k}"),
    ("expected overlap K x n / N", "141.9", f"{K * n / N:.1f}"),
    ("fold enrichment", "2.63", f"{k / (K * n / N):.2f}"),
    ("hypergeometric p", "1.95 x 10^-84", f"{10 ** (log10_p - exponent):.2f} x 10^{exponent}"),
], columns=["quantity", "Methods", "this table"])
checks["agrees"] = checks["Methods"] == checks["this table"]
with open(args.overlap) as f:
    step4 = json.load(f)
checks.loc[len(checks)] = [
    "N, K, n, k vs step 4 (overlap/methods_overlap.json)",
    f"{step4['N']:,}, {step4['K']:,}, {step4['n']:,}, {step4['k']}",
    f"{N:,}, {K:,}, {n:,}, {k}",
    (step4["N"], step4["K"], step4["n"], step4["k"]) == (N, K, n, k),
]

# ---- README ---------------------------------------------------------------------------------------------------------
built = datetime.date.today().isoformat()
tested_in = pair_table["expr_1_pct"].value_counts()
readme = pd.DataFrame([
    ("Table", "Pseudotime-associated genes (mutual information with pseudotime) in human and mouse fibroblasts and "
              "decidual cells, for strict one-to-one human-mouse orthologs. New version of human_mouse_mi (2 Jul 2026), "
              "rebuilt on the final pseudotime run."),
    ("Built", f"{built} by xmens_mi step 5 (scripts/supp_table.py) from the inputs below (with their sha256)."),
    ("Pseudotime", "obs['pseudotime'] of the final objects (rank-normalized diffusion pseudotime). Mouse: harmonized "
                   "pseudotime of 22 Oct 2026 (42,290 nuclei; the 907 apoptotic decidual nuclei have none). Human: 30 Oct "
                   "2026 refit on the fibroblast/decidual-only stromal scVI latent (47,533 nuclei)."),
    ("MI", "Per gene, MI between log-normalized expression and pseudotime (sklearn mutual_info_regression, n_neighbors = "
           "3, fixed random state); p from 1,000 pseudotime permutations, (count >= observed + 1) / 1,001; BH-adjusted "
           "q. MI is in bits."),
    ("Rows", f"Ensembl Compara 116 one-to-one pairs ({n_pairs:,}) with at least one gene tested: {len(pair_table):,} "
             f"pairs (tested in both {tested_in.get('both', 0):,}, human only {tested_in.get('human', 0):,}, mouse only "
             f"{tested_in.get('mouse', 0):,}); {n_neither_tested:,} pairs with neither gene tested are omitted. Sorted by "
             "human rank, then mouse rank."),
    ("Untested genes", "A gene detected in < 1% of nuclei was not tested and has no MI, p, q, rank or percentile."),
    ("human_gene / mouse_gene", "Ensembl Compara 116 one-to-one ortholog pair (homology type ortholog_one2one)."),
    ("top_10_pct_mi", "Species in which the gene is a top pseudotime-associated gene: MI above the 90th percentile of MI "
                      f"among that species' significant genes (mouse {mouse['threshold'] / LN2:.3f} bits; human "
                      f"{human['threshold'] / LN2:.3f} bits). both / human / mouse / neither. Genome-wide sets (1,426 "
                      "mouse, 1,489 human); within the background they give K, n and the overlap k."),
    ("human_mi_rank / mouse_mi_rank", "Rank of the gene's MI among all significant genes of its species, genome-wide "
                                      "(1 = highest). Empty if not significant."),
    ("human_mi_percentile / mouse_mi_percentile", "100 x (1 - (rank - 1) / number of significant genes in the "
                                                  "species). A top gene has percentile >= 90; the gene exactly at the "
                                                  "90th-percentile MI is not counted as top (strict >, as in the "
                                                  "Methods)."),
    ("human_mi / mouse_mi", "Mutual information with pseudotime (bits)."),
    ("human_q / mouse_q", "Benjamini-Hochberg q over all tested genes of the species."),
    ("human_p / mouse_p", "Permutation p-value (1,000 permutations; minimum 1/1,001)."),
    ("significant", "Species in which the gene is tested and q < 0.05: both / human / mouse / neither."),
    ("expr_1_pct", "Species in which the gene was detected in >= 1% of nuclei and so tested: both / human / mouse."),
    ("ortholog_univ", "TRUE if both genes are significant: the background of the cross-species comparison "
                      f"(N = {N:,} pairs; Methods hypergeometric test)."),
    ("Changes from the 2 Jul version", "Final pseudotime and MI of record; ortholog pairs from Ensembl Compara 116 (the "
                                       "July table had 15,134 pairs from the earlier homology set); untested genes now "
                                       "have no MI (the July run scored every gene and flagged expression afterwards); "
                                       "'significant' is now q < 0.05 among tested genes (the July file repeated "
                                       "expr_1_pct in this column); top_10_pct_mi uses the Methods threshold (MI above "
                                       "the 90th percentile among significant genes)."),
    ("Input: MI, mouse", f"{args.mi_mouse} (sha256 {sha256(args.mi_mouse)})"),
    ("Input: MI, human", f"{args.mi_human} (sha256 {sha256(args.mi_human)})"),
    ("Input: ortholog pairs", f"{args.pairs} (sha256 {sha256(args.pairs)})"),
    ("Sheet: overlap", "Front sheet. Live Excel formulas: change the cutoff in B4 and the thresholds, top-gene counts, "
                       "background N, K, n, overlap k, expected overlap, fold enrichment and the hypergeometric p "
                       "recalculate from the other sheets. At 90 it reproduces the Methods."),
    ("Sheet: human_genes / mouse_genes", "Every tested gene of the species (genome-wide, not only orthologs): MI (bits), "
                                         "p, q, significant (q < 0.05), rank among significant genes, and "
                                         "mi_if_significant (the column the percentile threshold is taken over)."),
    ("Sheet: p_value_terms", "Helper for the overlap sheet: log hypergeometric probability of each overlap x (log-gamma "
                             "form), so the upper tail P(X >= k) is summed without underflow (p near 1e-84 cannot be "
                             "computed as 1 - CDF)."),
    ("Sheet: checks", "The Methods numbers recomputed in Python when the workbook was built (static)."),
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
    checks.to_excel(writer, sheet_name="checks", index=False)
    wb = writer.book

    # Column widths, bold headers, frozen header row, and number formats of the pair table.
    widths = {
        "README": {"A": 34, "B": 140},
        "checks": {"A": 52, "B": 22, "C": 22, "D": 9},
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

    # Sheet overlap, the live front sheet: each row is (row, quantity, formula, Methods value, agreement test, format).
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
    overlap["A2"] = "Change the cutoff in B4; every value below recalculates. At 90 the values reproduce the Methods (column C)."
    rows = [
        (4, "Top-gene cutoff: percentile of MI among each species' significant genes", 90, None, None, "0"),
        (6, "Mouse MI threshold (bits)", f"=PERCENTILE({sig_m},$B$4/100)", 0.102, "=ROUND(B6,3)=C6", "0.000"),
        (7, "Human MI threshold (bits)", f"=PERCENTILE({sig_h},$B$4/100)", 0.039, "=ROUND(B7,3)=C7", "0.000"),
        (10, "Mouse genes tested", f"=COUNTA({genes_m})", 14261, "=B10=C10", "#,##0"),
        (11, "Human genes tested", f"=COUNTA({genes_h})", 15008, "=B11=C11", "#,##0"),
        (12, "Mouse significant genes (q < 0.05)", f"=COUNT({sig_m})", 14261, "=B12=C12", "#,##0"),
        (13, "Human significant genes (q < 0.05)", f"=COUNT({sig_h})", 14887, "=B13=C13", "#,##0"),
        (14, "Mouse top genes, genome-wide (MI > threshold)", f"=SUMPRODUCT(({sig_m}>B6)*1)", 1426, "=B14=C14", "#,##0"),
        (15, "Human top genes, genome-wide (MI > threshold)", f"=SUMPRODUCT(({sig_h}>B7)*1)", 1489, "=B15=C15", "#,##0"),
        (17, "Background pairs N (both genes significant)", f"=SUMPRODUCT(({univ}=TRUE)*1)", 9587, "=B17=C17", "#,##0"),
        (18, "Mouse top genes in the background, K", f"=SUMPRODUCT(({univ}=TRUE)*({mi_m}>B6))", 1133, "=B18=C18",
         "#,##0"),
        (19, "Human top genes in the background, n", f"=SUMPRODUCT(({univ}=TRUE)*({mi_h}>B7))", 1201, "=B19=C19",
         "#,##0"),
        (20, "Overlap k (top in both species)", f"=SUMPRODUCT(({univ}=TRUE)*({mi_m}>B6)*({mi_h}>B7))", 373,
         "=B20=C20", "#,##0"),
        (21, "Expected overlap, K x n / N", "=B18*B19/B17", 141.9, "=ROUND(B21,1)=C21", "0.0"),
        (22, "Expected overlap, % of the human top set", "=B21/B19", 0.12, "=ROUND(B22,2)=C22", "0%"),
        (23, "Observed overlap, % of the human top set", "=B20/B19", 0.31, "=ROUND(B23,2)=C23", "0%"),
        (24, "Fold enrichment, k / expected", "=B20/B21", 2.63, "=ROUND(B24,2)=C24", "0.00"),
        (25, "log10 p, hypergeometric P(X >= k)", "=p_value_terms!$F$2", None, None, "0.0000"),
        (26, "Hypergeometric p, P(X >= k): full-precision number (reads 0 below about 1E-308)", "=10^B25", None, None,
         "0.00E+00"),
        (27, "Hypergeometric p, as written in the Methods (text)", '=TEXT(10^(B25-INT(B25)),"0.00")&" x 10^"&INT(B25)',
         "1.95 x 10^-84", "=B27=C27", "@"),
    ]
    overlap["A3"], overlap["B3"], overlap["C3"], overlap["D3"] = "quantity", "value", "Methods (cutoff 90)", "agrees"
    for row, quantity, formula, methods_value, agreement, number_format in rows:
        overlap[f"A{row}"] = quantity
        overlap[f"B{row}"] = formula
        overlap[f"B{row}"].number_format = number_format
        if methods_value is not None:
            overlap[f"C{row}"] = methods_value
            overlap[f"C{row}"].number_format = number_format
        if agreement is not None:
            overlap[f"D{row}"] = f'=IF($B$4=90,{agreement[1:]},"(cutoff not 90)")'
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
    for column, width in {"A": 84, "B": 18, "C": 20, "D": 16}.items():
        overlap.column_dimensions[column].width = width

    # Sheet order, the front sheet first; Excel recalculates every formula when the file is opened.
    wb._sheets = [wb[sheet] for sheet in
                  ("overlap", "README", "human_mouse_mi", "human_genes", "mouse_genes", "p_value_terms", "checks")]
    wb.active = 0
    wb.calculation = CalcProperties(fullCalcOnLoad=True)

# ---- the CSVs and the provenance record; then the folder is renamed when complete -------------------------------------
pair_table.to_csv(os.path.join(partial, f"{NAME}__human_mouse_mi.csv"), index=False)
for species in SPECIES:
    gene_tables[species].to_csv(os.path.join(partial, f"{NAME}__{species}_genes.csv"), index=False)
checks.to_csv(os.path.join(partial, f"{NAME}_methods_checks.csv"), index=False)
provenance = {
    "built": built,
    "inputs": {key: {"path": path, "sha256": sha256(path)}
               for key, path in (("mi_mouse", args.mi_mouse), ("mi_human", args.mi_human), ("pairs", args.pairs))},
    "rows": len(pair_table),
    "omitted_neither_tested": n_neither_tested,
    "N": N, "K": K, "n": n, "k": k,
    "log10_p": log10_p,
    "all_checks_agree": bool(checks["agrees"].all()),
}
with open(os.path.join(partial, f"{NAME}_provenance.json"), "w") as f:
    json.dump(provenance, f, indent=1)
os.rename(partial, out_dir)
print(checks.to_string(index=False))
print(f"{len(pair_table):,} pairs ({n_neither_tested:,} with neither gene tested left out); wrote {out_dir}", flush=True)
