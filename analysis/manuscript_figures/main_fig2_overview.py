"""Main Figure 2: an exon-resolved, globally comparable map of disorder across human protein isoforms and protein classes.

a, toy gene with three isoforms: constitutive and alternative exons, shared and isoform-specific coding positions
(keyed by genomic codon) and the four IDR outcomes; b, exon view versus codon view in three splicing cases (cassette exon,
shifted splice site, different untranslated end); c, %IDR of canonical and alternative isoforms per protein class;
d, range of %IDR across each gene's isoforms per class; e, longest IDR vs %IDR for all isoforms with an IDR.
"""
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Patch
from scipy import stats
from style import *

CC = json.load(open(AN / 'manuscript_figures/outputs/cohort_counts.json')); RF = json.load(open(AN / 'reference_free/outputs/key_numbers.json'))
PC = AN / 'protein_classes/outputs'; G = pd.read_csv(PC / 'gene_annotations.csv', index_col=0); K = json.load(open(PC / 'key_numbers.json')); CL = K['classes_by_median_pct_idr']
iso = pd.read_csv(REPO / 'data/annotated/IDRisoforms_df_geneage.csv', low_memory=False, usecols=['base_accession', 'is_canonical', 'pct_idr', 'max_idr_len', 'n_idr_segments'])
iso['pclass'] = iso.base_accession.map(G.pclass)
CE = json.load(open(AN / 'reference_free/outputs/codon_vs_exon.json'))
fig = plt.figure(figsize=(7.1, 8.8))
gt = fig.add_gridspec(1, 2, top=0.985, bottom=0.715, left=0.08, right=0.99, width_ratios=[1.3, 1], wspace=0.1); gb = fig.add_gridspec(2, 2, top=0.665, bottom=0.04, hspace=0.5, wspace=0.3, width_ratios=[1, 1.35])

# ---- a: toy gene
ax = fig.add_subplot(gt[0, 0]); ax.set_xlim(0, 10.6); ax.set_ylim(-4.4, 3.2); ax.axis('off'); tag(ax, 'a', -0.12, 1.0)
EX = {1: (0.4, 1.4, [(0.0, 1.0, 0)]), 2: (2.0, 3.6, [(0, 0.45, 1), (0.45, 1, 0)]), 3: (4.2, 5.2, [(0, 1, 1)]), 4: (5.8, 6.6, [(0, 1, 1)]), 5: (7.2, 9.4, [(0, 0.3, 1), (0.3, 1, 0)])}
ISO = [('isoform 1', [1, 2, 3, 5]), ('isoform 2', [1, 3, 5]), ('isoform 3', [1, 2, 4, 5])]
ax.plot([0.2, 9.6], [2.6, 2.6], color=MUTED, lw=0.8, zorder=0)   # intron backbone, behind the exons
for e, (a, b, _) in EX.items():
    alt = sum(e in ex for _, ex in ISO) < 3
    ax.add_patch(Rectangle((a, 2.35), b - a, 0.5, fc='white', ec=ORANGE if alt else INK, lw=1.0 if alt else 0.7, zorder=2))
    ax.text((a + b) / 2, 2.6, 'A' if alt else 'C', ha='center', va='center', fontsize=6, color=ORANGE if alt else INK, zorder=3)
ax.text(0.2, 3.05, 'gene (exons on GRCh38): C, constitutive; A, alternative', fontsize=5.8, color=MUTED)
for k, (name, ex) in enumerate(ISO):
    y = 1.4 - 0.85 * k
    ax.text(-0.1, y + 0.2, name, ha='right', va='center', fontsize=5.8)
    for e in ex:
        a, b, parts = EX[e]
        for p0, p1, dis in parts:
            ax.add_patch(Rectangle((a + p0 * (b - a), y), (p1 - p0) * (b - a), 0.4, fc=BLUE if dis else LIGHT, ec='white', lw=0.4))
    xs = [(EX[e][0] + EX[e][1]) / 2 for e in ex]
    for e0, e1 in zip(ex[:-1], ex[1:]):
        ax.plot([EX[e0][1], EX[e1][0]], [y + 0.2, y + 0.2], color=MUTED, lw=0.4, ls=':')
# core / variable track
y = -1.35
for e, (a, b, _) in EX.items():
    n = sum(e in ex for _, ex in ISO)
    col = INK if n == 3 else (ORANGE if n == 1 else '#f2a37f')
    ax.add_patch(Rectangle((a, y), b - a, 0.3, fc=col, ec='white', lw=0.4))
    ax.text((a + b) / 2, y - 0.3, {3: 'shared', 2: 'majority', 1: 'minority'}[n], ha='center', va='center', fontsize=5.2, color=MUTED)
ax.text(-0.1, y + 0.15, 'positions', ha='right', va='center', fontsize=5.8)
ax.legend([Patch(color=BLUE), Patch(color=LIGHT)], ['IDR residues', 'ordered residues'], frameon=False, fontsize=5.6, loc='upper left', bbox_to_anchor=(-0.02, -0.02), ncol=2, handlelength=1)
# fates
ax.text(0.2, -2.35, 'IDR outcomes: an IDR of one isoform followed into another by its codons', fontsize=5.8, color=MUTED)
FT = [('intact', [(0, 1)]), ('internal deletion', [(0, 0.35), (0.65, 1)]), ('truncated', [(0, 0.6)]), ('removed', [])]
yy = -3.1
ax.text(-0.1, yy + 0.14, 'IDR in isoform X', ha='right', va='center', fontsize=5.4)
ax.text(-0.1, yy - 0.28, 'same codons, isoform Y', ha='right', va='center', fontsize=5.4)
for i, (f, keep) in enumerate(FT):
    x0 = 0.3 + i * 2.6
    ax.add_patch(Rectangle((x0, yy), 1.8, 0.28, fc=ORANGE, lw=0))
    ax.add_patch(Rectangle((x0, yy - 0.42), 1.8, 0.28, fc='none', ec=ORANGE, lw=0.5, ls=(0, (2, 1.5))))   # codons lost: dashed outline
    for k0, k1 in keep: ax.add_patch(Rectangle((x0 + k0 * 1.8, yy - 0.42), (k1 - k0) * 1.8, 0.28, fc=ORANGE, lw=0))
    ax.text(x0 + 0.9, yy - 0.85, f, ha='center', fontsize=5.6)

# ---- b: exon view vs codon view, three cases
ax = fig.add_subplot(gt[0, 1]); ax.set_xlim(0, 10); ax.set_ylim(-8.9, 0.5); ax.axis('off'); tag(ax, 'b', -0.02, 1.0)
X0 = 2.7   # track origin; row labels sit to the left
def ex(x0, x1, y, h=0.34, **k): ax.add_patch(Rectangle((X0 + x0, y - h / 2), x1 - x0, h, **k))
def gene(y, exons, utr=None):
    ax.plot([X0 + exons[0][1], X0 + exons[-1][0]], [y, y], color=MUTED, lw=0.5, zorder=0)
    for a, b in exons: ex(a, b, y, fc=LIGHT, ec=MUTED, lw=0.4)
    if utr: ex(utr[0], utr[1], y, h=0.16, fc=LIGHT, ec=MUTED, lw=0.4)
def view(y, segs):
    for a, b, spec in segs: ex(a, b, y, h=0.26, fc=ORANGE if spec else INK, lw=0)
E1, E3 = (0, 1.4), (5.6, 7.2)
CASES = [('1  cassette exon: both views agree', [E1, (2.6, 4.2), E3], None, [E1, E3], None,
          [(0, 1.4, 0), (2.6, 4.2, 1), (5.6, 7.2, 0)], [(0, 1.4, 0), (2.6, 4.2, 1), (5.6, 7.2, 0)]),
         ('2  shifted splice site: only the extra piece differs', [E1, (2.6, 4.8), E3], None, [E1, (2.6, 3.7), E3], None,
          [(0, 1.4, 0), (2.6, 4.8, 1), (5.6, 7.2, 0)], [(0, 1.4, 0), (2.6, 3.7, 0), (3.7, 4.8, 1), (5.6, 7.2, 0)]),
         ('3  different untranslated end: same protein', [E1, (2.6, 4.2), (5.6, 6.5)], (6.5, 7.3), [E1, (2.6, 4.2), (5.6, 6.5)], (6.5, 6.9),
          [(0, 1.4, 0), (2.6, 4.2, 0), (5.6, 6.5, 1)], [(0, 1.4, 0), (2.6, 4.2, 0), (5.6, 6.5, 0)])]
for i, (title, ea, ua, eb, ub, sv, cv) in enumerate(CASES):
    y0 = -2.95 * i
    ax.text(0, y0, title, fontsize=5.8, va='center')
    for dy, lab in [(-0.55, 'isoform A'), (-1.0, 'isoform B'), (-1.55, 'exon view'), (-1.95, 'codon view')]:
        ax.text(X0 - 0.25, y0 + dy, lab, fontsize=5.2, ha='right', va='center', color=INK if 'isoform' in lab else MUTED)
    gene(y0 - 0.55, ea, ua); gene(y0 - 1.0, eb, ub); view(y0 - 1.55, sv); view(y0 - 1.95, cv)
ax.text(X0 + 6.8, -6.95 + 0.55 + 0.2, 'UTR', fontsize=4.8, color=MUTED, ha='center', va='bottom')
ax.legend([Patch(color=INK), Patch(color=ORANGE)], ['constitutive exon / shared codon', 'alternative exon / isoform-specific codon'],
          frameon=False, fontsize=5.4, loc='lower left', bbox_to_anchor=(0.0, -0.03), ncol=1, handlelength=1)

# ---- c: violins per class
ax = fig.add_subplot(gb[0, :])
for i, c in enumerate(CL):
    for side, can, col in [(-1, True, BLUE), (1, False, ORANGE)]:
        v = iso[(iso.pclass == c) & (iso.is_canonical == can)].pct_idr.values
        kde = stats.gaussian_kde(v, bw_method=0.25); yv = np.linspace(0, 100, 200); dens = kde(yv); scale = 0.4 / dens.max(); d = scale * dens
        ax.fill_betweenx(yv, i, i + side * d, color=col, alpha=0.75, lw=0)
        q = np.percentile(v, [25, 50, 75]); half = scale * kde(q)       # violin half-width at each quantile: lines stop at the edge
        ax.plot([i, i + side * half[1]], [q[1], q[1]], color=INK, lw=1)
        for k in (0, 2): ax.plot([i, i + side * half[k]], [q[k], q[k]], color=INK, lw=0.5, ls=':')
ax.set_xticks(range(len(CL))); ax.set_xticklabels([f'{SHORT[c]}\n({int((G.pclass == c).sum()):,})' for c in CL], fontsize=6, rotation=30, ha='right')
ax.set_ylabel('residues in IDRs (%)'); ax.set_ylim(0, 100); ax.set_xlim(-0.6, len(CL) - 0.4); grid(ax)
ax.legend([Patch(color=BLUE), Patch(color=ORANGE)], ['canonical isoform', 'alternative isoforms'], frameon=False, fontsize=6.3, loc='upper left', ncol=2)
tag(ax, 'c', -0.06)
# ---- d: range across isoforms per class
ax = fig.add_subplot(gb[1, 0])
y = np.arange(len(CL))
for yi, c in zip(y, CL):
    v = G[(G.pclass == c)].idr_range.dropna().values; q = np.percentile(v, [10, 25, 50, 75, 90])
    ax.plot([q[0], q[4]], [yi, yi], color=MUTED, lw=0.7); ax.add_patch(plt.Rectangle((q[1], yi - 0.3), q[3] - q[1], 0.6, color=AQUA, alpha=0.75, lw=0))
    ax.plot([q[2]] * 2, [yi - 0.3, yi + 0.3], color=INK, lw=1)
ax.set_yticks(y); ax.set_yticklabels([SHORT[c] for c in CL], fontsize=6); ax.set_xlim(0, 60)
ax.set_xlabel("range of %IDR across a gene's isoforms\n(max − min, percentage points)"); grid(ax, 'x'); tag(ax, 'd', -0.75)
# ---- e: longest IDR vs %IDR
ax = fig.add_subplot(gb[1, 1])
w = iso[iso.n_idr_segments > 0]
ax.hexbin(np.log10(w.max_idr_len), w.pct_idr, gridsize=34, cmap='Blues', bins='log', mincnt=1, lw=0.1)
r = stats.spearmanr(w.max_idr_len, w.pct_idr)[0]
ax.text(0.03, 0.95, f'Spearman ρ = {r:.2f}\n{len(w):,} isoforms with an IDR', transform=ax.transAxes, va='top', fontsize=6)
ax.set_xticks(np.log10([20, 100, 1000])); ax.set_xticklabels(['20', '100', '1,000']); ax.set_xlabel('longest IDR (residues)'); ax.set_ylabel('% of residues that are disordered'); tag(ax, 'e', -0.15)
save(fig, 'main/fig2_overview')
