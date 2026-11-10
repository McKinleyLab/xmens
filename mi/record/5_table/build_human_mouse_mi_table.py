#!/usr/bin/env python
"""Build the human-mouse pseudotime mutual-information table (new version of results/tables/human_mouse_mi.csv) from the
MI of record of the final pseudotime run.

Inputs (cluster paths of record, 2 Nov 2026)
  --mi-mouse   xmens_final/results/mi_dpt/pt_harmonized_mouse_nall/mi_all.csv   (gene, mi [nats], p, q; tested genes only; written out in bits)
  --mi-human   xmens_final/results/mi_dpt/pt_harmonized_human_nall/mi_all.csv   (rerun 1 Nov 2026 on the 30 Oct human refit)
  --pairs      xmens_final/results/mi_dpt/overlap_of_record/ortholog_pairs_one2one.csv   (Ensembl Compara 116 one-to-one)
  --record-overlap  xmens_final/results/mi_dpt/overlap_of_record/methods_overlap.json      (optional, for the checks)
Definitions (Methods, "Pseudotime-associated genes and cross-species conservation")
  tested       gene detected in >= 1% of fibroblast/decidual nuclei (present in mi_all.csv)
  significant  tested and BH q < 0.05
  rank         rank of MI among the species' significant genes, genome-wide (1 = highest MI); percentile = 100 * (1 - (rank - 1) / N_significant)
  top          MI above the 90th percentile of MI among the species' significant genes (numpy percentile, linear interpolation)
  background   one-to-one pairs with both genes significant (ortholog_univ)
Outputs: <out>.xlsx (README, human_mouse_mi, checks), <out>.csv (data sheet), <out>_checks.csv, <out>_provenance.json
usage: build_human_mouse_mi_table.py --mi-mouse CSV --mi-human CSV --pairs CSV [--record-overlap JSON] --out PREFIX
"""
import argparse, datetime, hashlib, json
import numpy as np, pandas as pd
from scipy.stats import hypergeom

ap = argparse.ArgumentParser()
for name in ("--mi-mouse", "--mi-human", "--pairs", "--out"):
    ap.add_argument(name, required=True)
ap.add_argument("--record-overlap")
for name in ("--mi-mouse-label", "--mi-human-label", "--pairs-label"):
    ap.add_argument(name, help="path of record to show in README/provenance (default: the path given)")
a = ap.parse_args()
LABEL = {"mi_mouse": a.mi_mouse_label or a.mi_mouse, "mi_human": a.mi_human_label or a.mi_human, "pairs": a.pairs_label or a.pairs}
SPECIES = ("human", "mouse")
LN2 = float(np.log(2))
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

mi = {"human": pd.read_csv(a.mi_human), "mouse": pd.read_csv(a.mi_mouse)}
for s in SPECIES:
    assert not mi[s].gene.duplicated().any(), f"duplicate genes in {s} MI table"
    mi[s] = mi[s].set_index("gene")
pairs = pd.read_csv(a.pairs)[["human", "mouse"]]
n_pairs_input = len(pairs)
assert not pairs.human.duplicated().any() and not pairs.mouse.duplicated().any(), "pairs are not one-to-one"

stats = {}
for s in SPECIES:
    t = mi[s]
    sig = t[t.q < 0.05]
    thr = float(np.percentile(sig.mi, 90))
    rank = sig.mi.rank(ascending=False, method="first")
    stats[s] = dict(table=t, sig=sig, thr=thr, rank=rank, pct=100 * (1 - (rank - 1) / len(sig)),
                    top=set(sig.index[sig.mi > thr]))

T = pairs.rename(columns={"human": "human_gene", "mouse": "mouse_gene"}).copy()
flags = {}
for s in SPECIES:
    g = T[f"{s}_gene"]; st = stats[s]
    flags[s] = dict(tested=g.isin(st["table"].index).to_numpy(), sig=g.isin(st["sig"].index).to_numpy(),
                    top=g.isin(st["top"]).to_numpy())
    T[f"{s}_mi_rank"] = st["rank"].reindex(g).to_numpy()
    T[f"{s}_mi_percentile"] = st["pct"].reindex(g).to_numpy()
    T[f"{s}_mi"] = st["table"]["mi"].reindex(g).to_numpy() / LN2          # bits (the MI of record is in nats)
    for col in ("q", "p"):
        T[f"{s}_{col}"] = st["table"][col].reindex(g).to_numpy()
lab = lambda h, m: np.where(h & m, "both", np.where(h, "human", np.where(m, "mouse", "neither")))
T["top_10_pct_mi"] = lab(flags["human"]["top"], flags["mouse"]["top"])
T["significant"] = lab(flags["human"]["sig"], flags["mouse"]["sig"])
T["expr_1_pct"] = lab(flags["human"]["tested"], flags["mouse"]["tested"])
T["ortholog_univ"] = flags["human"]["sig"] & flags["mouse"]["sig"]
n_neither_tested = int((T.expr_1_pct == "neither").sum())
T = T[T.expr_1_pct != "neither"]
cols = ["human_gene", "mouse_gene", "top_10_pct_mi", "human_mi_rank", "mouse_mi_rank", "human_mi_percentile", "mouse_mi_percentile",
        "human_mi", "mouse_mi", "human_q", "mouse_q", "human_p", "mouse_p", "significant", "expr_1_pct", "ortholog_univ"]
T = T.sort_values(["human_mi_rank", "mouse_mi_rank", "human_gene"], na_position="last")[cols].reset_index(drop=True)
for s in SPECIES:
    T[f"{s}_mi_rank"] = T[f"{s}_mi_rank"].astype("Int64")

# ---- checks against the numbers stated in the Methods -----------------------------------------------------------------
U = T[T.ortholog_univ]
K = int(U.top_10_pct_mi.isin(["mouse", "both"]).sum()); n = int(U.top_10_pct_mi.isin(["human", "both"]).sum())
k = int((U.top_10_pct_mi == "both").sum()); N = len(U)
lp = float(hypergeom.logsf(k - 1, N, K, n) / np.log(10)); e10 = int(np.floor(lp))
checks = [("one-to-one pairs (Ensembl Compara 116)", "16,275", f"{n_pairs_input:,}"),
          ("genes tested, mouse / human", "14,261 / 15,008", f"{len(mi['mouse']):,} / {len(mi['human']):,}"),
          ("significant (q < 0.05), mouse / human", "14,261 / 14,887", f"{len(stats['mouse']['sig']):,} / {len(stats['human']['sig']):,}"),
          ("top genes genome-wide, mouse / human", "1,426 / 1,489", f"{len(stats['mouse']['top']):,} / {len(stats['human']['top']):,}"),
          ("top-gene MI threshold (bits), mouse / human", "0.102 / 0.039",
           f"{stats['mouse']['thr'] / np.log(2):.3f} / {stats['human']['thr'] / np.log(2):.3f}"),
          ("background pairs N (ortholog_univ = TRUE)", "9,587", f"{N:,}"),
          ("top mouse genes in background, K", "1,133", f"{K:,}"),
          ("top human genes in background, n", "1,201", f"{n:,}"),
          ("overlap k (top_10_pct_mi = both, in background)", "373", f"{k}"),
          ("expected overlap K x n / N", "141.9", f"{K * n / N:.1f}"),
          ("fold enrichment", "2.63", f"{k / (K * n / N):.2f}"),
          ("hypergeometric p", "1.95 x 10^-84", f"{10 ** (lp - e10):.2f} x 10^{e10}")]
C = pd.DataFrame(checks, columns=["quantity", "Methods", "this table"])
C["agrees"] = C["Methods"] == C["this table"]
if a.record_overlap:
    rec = json.load(open(a.record_overlap))
    C.loc[len(C)] = ["N, K, n, k vs pipeline overlap_of_record/methods_overlap.json",
                     f"{rec['N']:,}, {rec['K']:,}, {rec['n']:,}, {rec['k']}", f"{N:,}, {K:,}, {n:,}, {k}",
                     (rec["N"], rec["K"], rec["n"], rec["k"]) == (N, K, n, k)]

# ---- README sheet -------------------------------------------------------------------------------------------------------
built = datetime.date.today().isoformat()
counts = T.expr_1_pct.value_counts()
readme = [
 ("Table", "Pseudotime-associated genes (mutual information with pseudotime) in human and mouse fibroblasts and decidual cells, "
           "for strict one-to-one human-mouse orthologs. New version of human_mouse_mi (2 Jul 2026), rebuilt on the final pseudotime run."),
 ("Built", f"{built} by build_human_mouse_mi_table.py from the MI of record (inputs and sha256 below)."),
 ("Pseudotime", "obs['pseudotime'] of the final objects (rank-normalized diffusion pseudotime). Mouse: harmonized pseudotime of 22 Oct 2026 "
                "(42,290 nuclei; the 907 apoptotic decidual nuclei have none). Human: 30 Oct 2026 refit on the fibroblast/decidual-only "
                "stromal scVI latent (47,533 nuclei)."),
 ("MI", "Per gene, MI between log-normalized expression and pseudotime (sklearn mutual_info_regression, n_neighbors = 3, fixed random "
        "state); p from 1,000 pseudotime permutations, (count >= observed + 1) / 1,001; BH-adjusted q. MI is in bits."),
 ("Rows", f"Ensembl Compara 116 one-to-one pairs ({n_pairs_input:,}) with at least one gene tested: {len(T):,} pairs "
          f"(tested in both {counts.get('both', 0):,}, human only {counts.get('human', 0):,}, mouse only {counts.get('mouse', 0):,}); "
          f"{n_neither_tested:,} pairs with neither gene tested are omitted. Sorted by human rank, then mouse rank."),
 ("Untested genes", "A gene detected in < 1% of nuclei was not tested and has no MI, p, q, rank or percentile."),
 ("human_gene / mouse_gene", "Ensembl Compara 116 one-to-one ortholog pair (homology type ortholog_one2one)."),
 ("top_10_pct_mi", "Species in which the gene is a top pseudotime-associated gene: MI above the 90th percentile of MI among that species' "
                   f"significant genes (mouse {stats['mouse']['thr'] / LN2:.3f} bits; "
                   f"human {stats['human']['thr'] / LN2:.3f} bits). both / human / mouse / neither. "
                   "Genome-wide sets (1,426 mouse, 1,489 human); within the background they give K, n and the overlap k."),
 ("human_mi_rank / mouse_mi_rank", "Rank of the gene's MI among all significant genes of its species, genome-wide (1 = highest). Empty if not significant."),
 ("human_mi_percentile / mouse_mi_percentile", "100 x (1 - (rank - 1) / number of significant genes in the species). "
                   "A top gene has percentile >= 90; the gene exactly at the 90th-percentile MI is not counted as top (strict >, as in the Methods)."),
 ("human_mi / mouse_mi", "Mutual information with pseudotime (bits)."),
 ("human_q / mouse_q", "Benjamini-Hochberg q over all tested genes of the species."),
 ("human_p / mouse_p", "Permutation p-value (1,000 permutations; minimum 1/1,001)."),
 ("significant", "Species in which the gene is tested and q < 0.05: both / human / mouse / neither."),
 ("expr_1_pct", "Species in which the gene was detected in >= 1% of nuclei and so tested: both / human / mouse."),
 ("ortholog_univ", "TRUE if both genes are significant: the background of the cross-species comparison "
                   f"(N = {N:,} pairs; Methods hypergeometric test)."),
 ("Changes from the 2 Jul version", "Final pseudotime and MI of record; ortholog pairs from Ensembl Compara 116 (the July table had 15,134 pairs "
                   "from the earlier homology set); untested genes now have no MI (the July run scored every gene and flagged expression "
                   "afterwards); 'significant' is now q < 0.05 among tested genes (the July file repeated expr_1_pct in this column); "
                   "top_10_pct_mi uses the Methods threshold (MI above the 90th percentile among significant genes)."),
]
for key, k_, path in (("MI, mouse", "mi_mouse", a.mi_mouse), ("MI, human", "mi_human", a.mi_human), ("ortholog pairs", "pairs", a.pairs)):
    readme.append((f"Input: {key}", f"{LABEL[k_]} (sha256 {sha(path)})"))
R = pd.DataFrame(readme, columns=["field", "description"])

# ---- genome-wide gene sheets (all tested genes; feed the live thresholds on the overlap sheet) ------------------------------
G = {}
for s in SPECIES:
    t = mi[s].reset_index().sort_values("mi", ascending=False).reset_index(drop=True)
    t["significant"] = t.q < 0.05
    t["rank_among_significant"] = stats[s]["rank"].reindex(t.gene).astype("Int64").to_numpy()
    t["mi"] = t.mi / LN2                                                     # bits
    t["mi_if_significant"] = t.mi.where(t.significant)
    G[s] = t[["gene", "mi", "p", "q", "significant", "rank_among_significant", "mi_if_significant"]]

readme += [
 ("Sheet: overlap", "Front sheet. Live Excel formulas: change the cutoff in B4 and the thresholds, top-gene counts, background N, K, n, "
                    "overlap k, expected overlap, fold enrichment and the hypergeometric p recalculate from the other sheets. At 90 it "
                    "reproduces the Methods."),
 ("Sheet: human_genes / mouse_genes", "Every tested gene of the species (genome-wide, not only orthologs): MI (bits), p, q, significant "
                    "(q < 0.05), rank among significant genes, and mi_if_significant (the column the percentile threshold is taken over)."),
 ("Sheet: p_value_terms", "Helper for the overlap sheet: log hypergeometric probability of each overlap x (log-gamma form), so the upper "
                    "tail P(X >= k) is summed without underflow (p near 1e-84 cannot be computed as 1 - CDF)."),
 ("Sheet: checks", "The Methods numbers recomputed in Python when the workbook was built (static)."),
]
R = pd.DataFrame(readme, columns=["field", "description"])

# ---- write --------------------------------------------------------------------------------------------------------------
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.workbook.properties import CalcProperties
T.to_csv(f"{a.out}.csv", index=False)
C.to_csv(f"{a.out}_checks.csv", index=False)
for s in SPECIES:
    G[s].to_csv(f"{a.out}__{s}_genes.csv", index=False)
last = {"pairs": len(T) + 1, "human": len(G["human"]) + 1, "mouse": len(G["mouse"]) + 1}
NMAX = int(T.ortholog_univ.sum()) + 50                    # x runs 0..NMAX, beyond any min(K, n) the background allows
with pd.ExcelWriter(f"{a.out}.xlsx", engine="openpyxl") as xw:
    R.to_excel(xw, sheet_name="README", index=False)
    T.to_excel(xw, sheet_name="human_mouse_mi", index=False)
    G["human"].to_excel(xw, sheet_name="human_genes", index=False)
    G["mouse"].to_excel(xw, sheet_name="mouse_genes", index=False)
    C.to_excel(xw, sheet_name="checks", index=False)
    wb = xw.book
    bold = Font(bold=True)
    for name, widths in (("README", {"A": 34, "B": 140}), ("checks", {"A": 52, "B": 22, "C": 22, "D": 9}),
                         ("human_mouse_mi", {c_: 15 for c_ in "ABCDEFGHIJKLMNOP"}),
                         ("human_genes", {c_: 16 for c_ in "ABCDEFG"}), ("mouse_genes", {c_: 16 for c_ in "ABCDEFG"})):
        ws = wb[name]
        for col, w in widths.items():
            ws.column_dimensions[col].width = w
        for cell in ws[1]:
            cell.font = bold
        ws.freeze_panes = "A2"
    for row in wb["README"].iter_rows(min_row=2):
        row[1].alignment = Alignment(wrap_text=True, vertical="top"); row[0].alignment = Alignment(vertical="top")
    ws = wb["human_mouse_mi"]
    for col, f in {"F": "0.00", "G": "0.00", "H": "0.0000", "I": "0.0000", "J": "0.00E+00", "K": "0.00E+00", "L": "0.00E+00", "M": "0.00E+00"}.items():
        for cell in ws[col][1:]:
            cell.number_format = f

    # p_value_terms: log pmf of the hypergeometric overlap for x = 0..NMAX, then a log-sum-exp of the upper tail
    pt = wb.create_sheet("p_value_terms")
    N_, K_, n_, k_ = "overlap!$B$17", "overlap!$B$18", "overlap!$B$19", "overlap!$B$20"
    pt.append(["x", "log P(X = x), x in the upper tail [k, min(K, n)]", "exp(log P - max)", None, "max log P", "log10 p = log10 P(X >= k)"])
    lnC = lambda a_, b_: f"(GAMMALN({a_}+1)-GAMMALN({b_}+1)-GAMMALN({a_}-{b_}+1))"
    for i in range(NMAX + 1):
        r = i + 2
        x = f"A{r}"
        lp_ = f"{lnC(K_, x)}+{lnC(f'({N_}-{K_})', f'({n_}-{x})')}-{lnC(N_, n_)}"
        pt.append([i, f"=IF(AND({x}>={k_},{x}<=MIN({K_},{n_}),{n_}-{x}<={N_}-{K_}),{lp_},-1E+300)", f"=EXP(B{r}-$E$2)"])
    pt["E2"] = f"=MAX($B$2:$B${NMAX + 2})"
    pt["F2"] = f"=($E$2+LN(SUM($C$2:$C${NMAX + 2})))/LN(10)"
    for cell in pt[1]:
        cell.font = bold
    pt.column_dimensions["B"].width = 44; pt.column_dimensions["C"].width = 18; pt.column_dimensions["E"].width = 14; pt.column_dimensions["F"].width = 26

    # overlap: the live front sheet
    ov = wb.create_sheet("overlap", 0)
    P_, H_, I_ = (f"human_mouse_mi!${c_}$2:${c_}${last['pairs']}" for c_ in ("P", "H", "I"))
    FH, FM = f"human_genes!$G$2:$G${last['human']}", f"mouse_genes!$G$2:$G${last['mouse']}"
    AH, AM = f"human_genes!$A$2:$A${last['human']}", f"mouse_genes!$A$2:$A${last['mouse']}"
    ov["A1"] = "Cross-species overlap of top pseudotime-associated genes (live: computed from the other sheets of this workbook)"
    ov["A2"] = "Change the cutoff in B4; every value below recalculates. At 90 the values reproduce the Methods (column C)."
    rows = [
        (4, "Top-gene cutoff: percentile of MI among each species' significant genes", 90, None, None, "0"),
        (6, "Mouse MI threshold (bits)", f"=PERCENTILE({FM},$B$4/100)", 0.102, "=ROUND(B6,3)=C6", "0.000"),
        (7, "Human MI threshold (bits)", f"=PERCENTILE({FH},$B$4/100)", 0.039, "=ROUND(B7,3)=C7", "0.000"),
        (10, "Mouse genes tested", f"=COUNTA({AM})", 14261, "=B10=C10", "#,##0"),
        (11, "Human genes tested", f"=COUNTA({AH})", 15008, "=B11=C11", "#,##0"),
        (12, "Mouse significant genes (q < 0.05)", f"=COUNT({FM})", 14261, "=B12=C12", "#,##0"),
        (13, "Human significant genes (q < 0.05)", f"=COUNT({FH})", 14887, "=B13=C13", "#,##0"),
        (14, "Mouse top genes, genome-wide (MI > threshold)", f"=SUMPRODUCT(({FM}>B6)*1)", 1426, "=B14=C14", "#,##0"),
        (15, "Human top genes, genome-wide (MI > threshold)", f"=SUMPRODUCT(({FH}>B7)*1)", 1489, "=B15=C15", "#,##0"),
        (17, "Background pairs N (both genes significant)", f"=SUMPRODUCT(({P_}=TRUE)*1)", 9587, "=B17=C17", "#,##0"),
        (18, "Mouse top genes in the background, K", f"=SUMPRODUCT(({P_}=TRUE)*({I_}>B6))", 1133, "=B18=C18", "#,##0"),
        (19, "Human top genes in the background, n", f"=SUMPRODUCT(({P_}=TRUE)*({H_}>B7))", 1201, "=B19=C19", "#,##0"),
        (20, "Overlap k (top in both species)", f"=SUMPRODUCT(({P_}=TRUE)*({I_}>B6)*({H_}>B7))", 373, "=B20=C20", "#,##0"),
        (21, "Expected overlap, K x n / N", "=B18*B19/B17", 141.9, "=ROUND(B21,1)=C21", "0.0"),
        (22, "Expected overlap, % of the human top set", "=B21/B19", 0.12, "=ROUND(B22,2)=C22", "0%"),
        (23, "Observed overlap, % of the human top set", "=B20/B19", 0.31, "=ROUND(B23,2)=C23", "0%"),
        (24, "Fold enrichment, k / expected", "=B20/B21", 2.63, "=ROUND(B24,2)=C24", "0.00"),
        (25, "log10 p, hypergeometric P(X >= k)", "=p_value_terms!$F$2", None, None, "0.0000"),
        (26, "Hypergeometric p, P(X >= k): full-precision number (reads 0 below about 1E-308)", "=10^B25", None, None, "0.00E+00"),
        (27, "Hypergeometric p, as written in the Methods (text)", '=TEXT(10^(B25-INT(B25)),"0.00")&" x 10^"&INT(B25)', "1.95 x 10^-84", "=B27=C27", "@"),
    ]
    ov["A3"], ov["B3"], ov["C3"], ov["D3"] = "quantity", "value", "Methods (cutoff 90)", "agrees"
    for r, label, formula, methods, agree, fmt in rows:
        ov[f"A{r}"] = label; ov[f"B{r}"] = formula; ov[f"B{r}"].number_format = fmt
        if methods is not None:
            ov[f"C{r}"] = methods; ov[f"C{r}"].number_format = fmt
        if agree is not None:
            ov[f"D{r}"] = f'=IF($B$4=90,{agree[1:]},"(cutoff not 90)")'
    ov["A29"] = ("N, K, n and k are counted in human_mouse_mi; the thresholds and genome-wide counts in human_genes / mouse_genes; "
                 "p is summed in p_value_terms. MI thresholds use PERCENTILE (linear interpolation, as numpy.percentile).")
    for cell in ov[3]:
        cell.font = bold
    ov["A1"].font = Font(bold=True, size=12)
    ov["B4"].fill = PatternFill("solid", fgColor="FFF2CC"); ov["B4"].font = bold
    for r in range(17, 28):                                   # the overlap-test results: N, K, n, k, expected, %, fold, p
        ov[f"B{r}"].font = bold
    ov.column_dimensions["A"].width = 84; ov.column_dimensions["B"].width = 18; ov.column_dimensions["C"].width = 20; ov.column_dimensions["D"].width = 16
    wb._sheets = [wb[n] for n in ("overlap", "README", "human_mouse_mi", "human_genes", "mouse_genes", "p_value_terms", "checks")]
    wb.active = 0
    wb.calculation = CalcProperties(fullCalcOnLoad=True)
json.dump(dict(built=built, inputs={k: dict(path=LABEL[k], sha256=sha(v)) for k, v in
                                    dict(mi_mouse=a.mi_mouse, mi_human=a.mi_human, pairs=a.pairs).items()},
               rows=len(T), omitted_neither_tested=n_neither_tested, N=N, K=K, n=n, k=k, log10_p=lp,
               all_checks_agree=bool(C.agrees.all())), open(f"{a.out}_provenance.json", "w"), indent=1)
print(C.to_string(index=False)); print("rows", len(T), "| omitted (neither tested)", n_neither_tested)
