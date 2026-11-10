"""Supplementary figure (draft): the two recurrence cutoffs behind the floor-6 consensus regulons.

  a  runs (of 100) in which each TF formed a regulon; bars at >= 80 runs (kept TFs) in species colour
  b  regulons kept vs TF threshold (with the >= 10-target size check at floor 6), per species and in both species
     (TFs paired by Ensembl one-to-one orthologs, as in the quadrant figure)
  c  share of TF-target pairs, and of their occurrence weight, kept vs target floor, within kept regulons
  d  regulon size (median and IQR) vs target floor; dotted = median single-run regulon size of the same TFs

Inputs (artifacts)
  report_occ6_{mouse,human}.tsv             ae8e54fe / 6e5a2f8e  per-TF report of the floor-6 consensus: n_present (runs
                                                                  with a regulon), n_targets_union (distinct targets over
                                                                  all runs), n_targets_kept (targets in >= 6 runs), kept
  {mouse,human}_consensus_regulons_long.tsv cae18a95 / 056da4ca  every pair of a kept TF recurring in >= 2 runs, with
                                                                  n_seeds; pairs seen in 1 run = n_targets_union - rows
  perrun_sizes.json                         c0970b1c             per-TF median regulon size within single runs
  ortholog_pairs_one2one.csv                f4b0ede2             Ensembl one-to-one human-mouse symbol pairs
Occurrence weight of a pair = runs containing the pair / runs in which its TF formed a regulon (the AUCell weight).
Output: regulon_recurrence_cutoffs.svg / .png (manuscript Fig. S28 layout; 2026-11-02: 'TF gate' -> 'TF threshold',
'target floor' -> 'target threshold' in the axis labels, per S. Eddy's comments)
"""
import json
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

REP6 = {'mouse': 'ae8e54fe-6eab-4ce1-a9b3-2cca300a0116', 'human': '6e5a2f8e-de95-4cd8-954a-5df5e2491988'}
LONG = {'mouse': 'cae18a95-9437-4665-b002-33690a40a4f9', 'human': '056da4ca-67ae-4cf0-aa01-108e0d75a626'}
PERRUN, ORTHO = 'c0970b1c-578d-4b39-a683-9f4ea2a8c83c', 'f4b0ede2-d511-4e32-a029-0360064243aa'
GATE, FLOOR, MIN_T = 80, 6, 10
BASE, SMALL = 16, 14
COL = {'mouse': '#4d4d4d', 'human': '#1b9e77'}
GREY, REF = '#c8c8c8', '0.55'
OUT = globals().get('OUT', 'regulon_recurrence_cutoffs')

# ---- data
R = {sp: pd.read_csv(host.artifact_path(v), sep='\t').set_index('name') for sp, v in REP6.items()}
L = {sp: pd.read_csv(host.artifact_path(v), sep='\t') for sp, v in LONG.items()}
PR = json.load(open(host.artifact_path(PERRUN)))
O = pd.read_csv(host.artifact_path(ORTHO))[['human', 'mouse']].drop_duplicates()
F = np.arange(1, 101)
D = {}
for sp in ('mouse', 'human'):
    r = R[sp]
    assert r.index.str.endswith('(+)').all()
    k = r[r.kept]
    assert len(k) == int(((r.n_present >= GATE) & (r.n_targets_kept >= MIN_T)).sum()) == int((r.n_present >= GATE).sum())
    d = L[sp]
    assert set(d.regulon) == set(k.index)
    single = k.n_targets_union - d.groupby('regulon').size().reindex(k.index)
    assert (single >= 0).all()
    assert (d[d.n_seeds >= FLOOR].groupby('regulon').size().reindex(k.index) == k.n_targets_kept).all()
    ns = d.n_seeds.to_numpy()
    wv = (d.n_seeds / d.regulon.map(k.n_present)).to_numpy()
    w1 = float((single / k.n_present).sum())
    tot_n, tot_w = len(d) + int(single.sum()), wv.sum() + w1
    D[sp] = dict(
        r=r, k=k, tot_n=tot_n,
        pairs=np.array([((ns >= f).sum() + (single.sum() if f == 1 else 0)) / tot_n for f in F]),
        weight=np.array([(wv[ns >= f].sum() + (w1 if f == 1 else 0)) / tot_w for f in F]),
        sizes=np.vstack([k.n_targets_union.to_numpy() if f == 1 else
                         d[d.n_seeds >= f].groupby('regulon').size().reindex(k.index).fillna(0).to_numpy() for f in F]),
        perrun=float(np.median([PR[sp]['per_tf_median_size'][n] for n in k.index])),
        gate=np.array([int(((r.n_present >= g) & (r.n_targets_kept >= MIN_T)).sum()) for g in F]))


def passing(sp, g):
    r = R[sp]
    return {n[:-3] for n in r.index[(r.n_present >= g) & (r.n_targets_kept >= MIN_T)]}


PAIRED = np.array([int((O.human.isin(passing('human', g)) & O.mouse.isin(passing('mouse', g))).sum()) for g in F])
assert PAIRED[GATE - 1] == 111, PAIRED[GATE - 1]

# ---- figure
apply_figure_style(sizes=(BASE, SMALL, SMALL))
_have = {f.name for f in mpl.font_manager.fontManager.ttflist}
mpl.rcParams.update({'font.family': [f for f in ('Helvetica Neue', 'Helvetica', 'Arial', 'Nimbus Sans', 'DejaVu Sans') if f in _have]})
plt.close('all')
fig = plt.figure(figsize=(11.0, 8.4), dpi=300)
outer = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.30, left=0.085, right=0.975, top=0.925, bottom=0.085)
ga = outer[0, 0].subgridspec(2, 1, hspace=0.22)
axa1 = fig.add_subplot(ga[0]); axa2 = fig.add_subplot(ga[1], sharex=axa1, sharey=axa1)
axb, axc, axd = fig.add_subplot(outer[0, 1]), fig.add_subplot(outer[1, 0]), fig.add_subplot(outer[1, 1])

# a: TF recurrence
ymax = max(np.bincount(D[sp]['r'].n_present, minlength=101)[1:].max() for sp in D) * 1.32
for ax, sp in ((axa1, 'mouse'), (axa2, 'human')):
    cnt = np.bincount(D[sp]['r'].n_present.to_numpy(), minlength=101)[1:]
    ax.bar(F, cnt, width=1.0, color=np.where(F >= GATE, COL[sp], GREY), linewidth=0)
    ax.axvline(GATE - 0.5, color=REF, lw=0.9)
    ax.set_ylim(0, ymax); ax.set_xlim(0, 101.5)
    ax.text(0.015, 0.95, f'{sp.capitalize()}: {int(cnt.sum())} TFs with a regulon\nin at least one run',
            transform=ax.transAxes, ha='left', va='top', fontsize=SMALL, color=COL[sp])
    ax.text(GATE - 2, ymax * 0.40, f'{len(D[sp]["k"])} kept', ha='right', va='center', fontsize=SMALL, color=COL[sp])
    ax.set_ylabel('TFs')
plt.setp(axa1.get_xticklabels(), visible=False)
axa2.set_xticks([1, 20, 40, 60, 80, 100])
axa2.set_xlabel('Runs (of 100) in which the TF formed a regulon')

# b: regulons kept vs gate
for sp in ('mouse', 'human'):
    axb.plot(F, D[sp]['gate'], color=COL[sp], lw=1.8, label=sp.capitalize())
axb.plot(F, PAIRED, color='black', lw=1.4, ls=(0, (4, 2)), label='In both species')
axb.axvline(GATE, color=REF, lw=0.9)
for v, c_, dx, dy, ha, va in ((D['mouse']['gate'][GATE - 1], COL['mouse'], 1.5, 6, 'left', 'bottom'),
                              (D['human']['gate'][GATE - 1], COL['human'], -1.5, 17, 'right', 'bottom'),
                              (PAIRED[GATE - 1], 'black', -1.5, -6, 'right', 'top')):
    axb.plot([GATE], [v], 'o', ms=4.5, color=c_, zorder=4)
    axb.text(GATE + dx, v + dy, f'{v}', ha=ha, va=va, fontsize=SMALL, color=c_)
axb.set_xlim(0, 101.5); axb.set_ylim(0, max(D['mouse']['gate'].max(), D['human']['gate'].max()) * 1.12)
axb.set_xticks([1, 20, 40, 60, 80, 100])
axb.set_xlabel('TF threshold (minimum runs, of 100)'); axb.set_ylabel('Regulons kept')
axb.legend(loc='lower left', frameon=False, fontsize=SMALL, handlelength=2.2)

# c: pairs and weight kept vs floor
LOGT = [1, 2, 5, 10, 20, 50, 100]
for sp in ('mouse', 'human'):
    axc.plot(F, 100 * D[sp]['pairs'], color=COL[sp], lw=1.8)
    axc.plot(F, 100 * D[sp]['weight'], color=COL[sp], lw=1.8, ls=(0, (4, 2)))
    for kind in ('pairs', 'weight'):
        y = 100 * D[sp][kind][FLOOR - 1]
        fx_, dy, ha, va = (1.1, 2.5, 'left', 'bottom') if sp == 'mouse' else (1 / 1.1, -2.0, 'right', 'top')
        axc.plot([FLOOR], [y], 'o', ms=4.5, color=COL[sp], zorder=4)
        axc.text(FLOOR * fx_, y + dy, f'{y:.0f} %', ha=ha, va=va, fontsize=SMALL, color=COL[sp])
axc.axvline(FLOOR, color=REF, lw=0.9)
axc.set_xscale('log'); axc.set_xticks(LOGT); axc.set_xticklabels([str(t) for t in LOGT])
axc.xaxis.set_minor_locator(mpl.ticker.NullLocator()); axc.set_xlim(0.92, 108); axc.set_ylim(0, 104)
axc.set_xlabel('Target threshold (minimum runs)'); axc.set_ylabel('Kept (%)')
axc.text(0.98, 0.97, '100 % = all pairs\nseen in any run', transform=axc.transAxes,
         ha='right', va='top', fontsize=SMALL, color='0.35')
axc.legend(handles=[mpl.lines.Line2D([], [], color='0.3', lw=1.8, label='TF\u2013target pairs'),
                    mpl.lines.Line2D([], [], color='0.3', lw=1.8, ls=(0, (4, 2)), label='Occurrence weight')],
           loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=SMALL, handlelength=2.2,
           columnspacing=1.6, borderaxespad=0.2)

# d: regulon size vs floor
for sp in ('mouse', 'human'):
    s = D[sp]['sizes']
    med = np.median(s, axis=1); q1, q3 = np.percentile(s, [25, 75], axis=1)
    ok = med > 0
    axd.fill_between(F[ok], np.clip(q1[ok], 1, None), q3[ok], color=COL[sp], alpha=0.16, lw=0)
    axd.plot(F[ok], med[ok], color=COL[sp], lw=1.8)
    axd.axhline(D[sp]['perrun'], color=COL[sp], lw=1.1, ls=':')
    m6 = med[FLOOR - 1]
    axd.plot([FLOOR], [m6], 'o', ms=4.5, color=COL[sp], zorder=4)
    if sp == 'mouse':
        axd.text(FLOOR * 1.12, m6 * 1.05, f'{m6:.0f}', ha='left', va='bottom', fontsize=SMALL, color=COL[sp])
    else:
        axd.text(FLOOR / 1.12, m6 / 1.05, f'{m6:.0f}', ha='right', va='top', fontsize=SMALL, color=COL[sp])
axd.text(1.0, D['mouse']['perrun'] * 1.08, f'one run: {D["mouse"]["perrun"]:.0f}', ha='left', va='bottom', fontsize=SMALL, color=COL['mouse'])
axd.text(1.0, D['human']['perrun'] / 1.08, f'one run: {D["human"]["perrun"]:.0f}', ha='left', va='top', fontsize=SMALL, color=COL['human'])
axd.axvline(FLOOR, color=REF, lw=0.9)
axd.set_xscale('log'); axd.set_yscale('log')
axd.set_xticks(LOGT); axd.set_xticklabels([str(t) for t in LOGT]); axd.xaxis.set_minor_locator(mpl.ticker.NullLocator())
axd.set_xlim(0.92, 108); axd.set_ylim(1, 6000)
axd.set_yticks([1, 10, 100, 1000]); axd.set_yticklabels(['1', '10', '100', '1k']); axd.yaxis.set_minor_locator(mpl.ticker.NullLocator())
axd.set_xlabel('Target threshold (minimum runs)'); axd.set_ylabel('Targets per regulon')

for ax in fig.axes:                                   # bold axis labels, as in the manuscript
    ax.xaxis.label.set_fontweight('bold'); ax.yaxis.label.set_fontweight('bold')
for ax, lt, dx in ((axa1, 'A', -0.13), (axb, 'B', -0.13), (axc, 'C', -0.13), (axd, 'D', -0.13)):
    panel_letter(ax, lt, dx=-0.15, dy=1.07, case='upper', fontsize=BASE + 4)

fig.canvas.draw()
Rr = fig.canvas.get_renderer()
fig.savefig(OUT + '.png')
with mpl.rc_context({'svg.fonttype': 'none'}):
    fig.savefig(OUT + '.svg')

# ---- audit (figure-style 9.1)
fig.canvas.draw()
OFFVIEW = set()   # tick labels at ticks outside the view limits are never drawn
for a in fig.axes:
    for axis, lim in ((a.xaxis, a.get_xlim()), (a.yaxis, a.get_ylim())):
        lo, hi = sorted(lim)
        for tk in axis.get_major_ticks() + axis.get_minor_ticks():
            loc = tk.get_loc()
            if loc is None or not (lo - 1e-9 <= loc <= hi + 1e-9):
                OFFVIEW.update({tk.label1, tk.label2})
_offview = lambda t: t in OFFVIEW
texts = [(t, t.get_window_extent(Rr)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible() and not _offview(t)]
spines = [(s, s.get_window_extent(Rr)) for a in fig.axes for s in a.spines.values() if s.get_visible()]
tick = {a: set(a.get_xticklabels(which='both') + a.get_yticklabels(which='both')) for a in fig.axes}
OVERLAPS = [(a.get_text()[:30], b.get_text()[:30]) for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]
OVERLAPS += [(t.get_text()[:30], 'spine') for t, bt in texts for s, bs in spines if bt.overlaps(bs) and t not in tick[s.axes]]
fb = fig.bbox
OUTSIDE = [t.get_text()[:30] for t, b in texts if b.x0 < fb.x0 or b.x1 > fb.x1 or b.y0 < fb.y0 or b.y1 > fb.y1]
