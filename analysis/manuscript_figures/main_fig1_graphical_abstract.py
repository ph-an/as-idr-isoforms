"""Main Figure 1: graphical abstract.

a, the question: one gene, exons that encode ordered (green) or disordered (pink) sequence, and the isoforms that alternative
splicing makes from them, placed on an order-disorder axis; b, the approach in five steps; c, the five main findings, each with
its headline numbers (drawn from the analysis outputs, not typed in).
"""
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch
from style import *

GREEN, GREEN_E, PINK, PINK_E = '#bfe0a8', '#6aa84f', '#f6b9b9', '#e06666'
CC = json.load(open(AN / 'manuscript_figures/outputs/cohort_counts.json')); RF = json.load(open(AN / 'reference_free/outputs/key_numbers.json'))
PC = json.load(open(AN / 'protein_classes/outputs/key_numbers.json')); EV = json.load(open(AN / 'evolution/outputs/key_numbers.json'))
EA = json.load(open(AN / 'exon_age/outputs/key_numbers.json'))
A = pd.read_csv(AN / 'protein_classes/outputs/A_pct_idr_by_class.csv').set_index('pclass'); C = pd.read_csv(AN / 'protein_classes/outputs/C_exon_coupling_by_class.csv').set_index('pclass')
CL = PC['classes_by_median_pct_idr']

fig = plt.figure(figsize=(7.1, 8.5))
gs = fig.add_gridspec(3, 1, height_ratios=[1.43, 0.42, 0.8], hspace=0.3, top=0.985, bottom=0.05, left=0.01, right=0.99)

def box(ax, x, y, w, h, fc, ec, lw=0.8, r=0.06, z=2, ls='-'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec=ec, lw=lw, zorder=z, ls=ls))

def arrow(ax, p, q, c=MUTED, lw=0.8, ls='-', ms=7):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=ms, color=c, lw=lw, ls=ls, zorder=1))

def wiggle(ax, x0, y0, length, amp=0.07, n=3.5, ang=0.0, c=PINK_E, lw=1.1, z=3):
    t = np.linspace(0, 1, 60); xs = t * length; ys = amp * np.sin(2 * np.pi * n * t)
    ca, sa = np.cos(ang), np.sin(ang); X = x0 + xs * ca - ys * sa; Y = y0 + xs * sa + ys * ca
    ax.plot(X, Y, color=c, lw=lw, solid_capstyle='round', zorder=z); return X[-1], Y[-1]

def protein(ax, x, y, parts, s=1.0):
    """parts: sequence of 'D' (folded domain) and 'I' (disordered stretch), drawn left to right."""
    cx = x
    for i, p in enumerate(parts):
        if p == 'D':
            r = 0.17 * s; ax.add_patch(Circle((cx + r, y), r, fc=GREEN, ec=GREEN_E, lw=0.8, zorder=3)); cx += 2 * r
        else:
            ang = 0.5 if i % 2 else -0.5; ex, ey = wiggle(ax, cx, y, 0.42 * s, ang=ang * 0.4); cx = ex
    return cx

# ---------------- a: concept
ax = fig.add_subplot(gs[0]); ax.set_xlim(0, 10); ax.set_ylim(-1.3, 6.2); ax.axis('off'); tag(ax, 'a', 0.0, 0.97)
ax.text(0.35, 5.95, 'Alternative splicing decides which ordered and disordered pieces of a gene each protein carries',
        fontsize=8, fontweight='bold', va='center')
# gene
EX = [('E1', 'o'), ('E2', 'd'), ('E3', 'o'), ('E4', 'd'), ('E5', 'd')]
gx = np.linspace(1.6, 8.6, 5); ax.plot([1.0, 9.2], [5.2, 5.2], color=MUTED, lw=1, zorder=0)
ax.text(0.35, 5.2, 'gene', fontsize=6.5, va='center', color=MUTED); ax.text(9.35, 5.2, "3′", fontsize=6, va='center', color=MUTED); ax.text(0.85, 5.2, "5′", fontsize=6, va='center', ha='right', color=MUTED)
for (lab, k), x in zip(EX, gx):
    box(ax, x - 0.45, 5.02, 0.9, 0.36, GREEN if k == 'o' else PINK, GREEN_E if k == 'o' else PINK_E, r=0.05)
    ax.text(x, 5.2, lab, ha='center', va='center', fontsize=6.2, zorder=4)
# unresolved exon callout (from the poster): an exon can encode both
ax.text(5.0, 4.62, 'alternative splicing', ha='center', fontsize=6.5, color=INK, fontweight='bold', bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='none'))
ISO = [['E1', 'E3'], ['E1', 'E2', 'E3'], ['E1', 'E2', 'E3', 'E5'], ['E1', 'E2', 'E4', 'E5'], ['E2', 'E4', 'E5']]
PROT = [['D', 'D'], ['D', 'I', 'D'], ['D', 'I', 'D', 'I'], ['D', 'I', 'I', 'I'], ['I', 'I', 'I']]
kind = dict(EX); xs = np.linspace(1.35, 8.85, 5)
for k, (ex, pr, x) in enumerate(zip(ISO, PROT, xs)):
    arrow(ax, (5.0, 4.45), (x + 0.35, 4.0), c=LIGHT, lw=0.6, ms=5)
    w = 0.34; x0 = x + 0.35 - w * len(ex) / 2
    for j, e in enumerate(ex):
        o = kind[e] == 'o'; box(ax, x0 + j * w, 3.62, w - 0.02, 0.3, GREEN if o else PINK, GREEN_E if o else PINK_E, lw=0.6, r=0.03)
        ax.text(x0 + j * w + (w - 0.02) / 2, 3.77, e, ha='center', va='center', fontsize=4.6, zorder=4)
    arrow(ax, (x + 0.35, 3.52), (x + 0.35, 3.15), c=LIGHT, lw=0.6, ms=5)
    protein(ax, x - 0.15, 2.72, pr, s=0.95)
ax.text(0.2, 3.77, 'mRNAs', fontsize=6.5, va='center', color=MUTED); ax.text(0.2, 2.72, 'proteins', fontsize=6.5, va='center', color=MUTED)
# order-disorder axis
from matplotlib.colors import LinearSegmentedColormap
cm = LinearSegmentedColormap.from_list('od', [GREEN, '#ffffff', PINK])
ax.imshow(np.linspace(0, 1, 256)[None, :], extent=(1.0, 9.3, 1.75, 1.97), cmap=cm, aspect='auto', zorder=1)
for t, xx in [('fully ordered', 1.0), ('mostly ordered', 3.08), ('50%', 5.15), ('mostly disordered', 7.23), ('fully disordered', 9.3)]:
    ax.text(xx, 1.6, t, ha='center', va='top', fontsize=5.6, color=MUTED)
for x, f in zip(xs, [0.0, 0.3, 0.45, 0.72, 1.0]):
    ax.plot([x + 0.35, 1.0 + 8.3 * f], [2.45, 2.0], color=LIGHT, lw=0.5, ls=':')
# legend for protein parts
box(ax, 0.4, 0.98, 0.35, 0.22, GREEN, GREEN_E, r=0.03); ax.text(0.85, 1.09, 'exon encoding ordered sequence', fontsize=5.8, va='center')
box(ax, 0.4, 0.62 - 0.11, 0.35, 0.22, PINK, PINK_E, r=0.03); ax.text(0.85, 0.62, 'exon encoding disordered sequence', fontsize=5.8, va='center')
protein(ax, 3.3, 1.09, ['D'], s=0.8); ax.text(3.75, 1.09, 'folded domain', fontsize=5.8, va='center')
wiggle(ax, 3.3, 0.62, 0.35); ax.text(3.75, 0.62, 'disordered region', fontsize=5.8, va='center')
# exon zoom (from the poster): an exon is a run of residues, each predicted ordered or disordered
ax.text(0.4, 0.14, 'one exon, residue by residue', fontsize=5.8, va='center', color=MUTED)
ZX0, ZW, ZY = 0.4, 4.4, -0.62; nres = 20; dis = np.array([0] * 6 + [1] * 11 + [0] * 3, bool)
box(ax, 2.25, -0.02, 0.7, 0.26, PINK, 'none', lw=0, r=0.03, z=2)   # mixed exon: pink with green stripes, as on the poster
plt.rcParams['hatch.linewidth'] = 2.2
ax.add_patch(FancyBboxPatch((2.25, -0.02), 0.7, 0.26, boxstyle='round,pad=0,rounding_size=0.03', fc='none', ec=GREEN, lw=0, hatch='///', zorder=3))
box(ax, 2.25, -0.02, 0.7, 0.26, 'none', DIS_E if 'DIS_E' in globals() else PINK_E, lw=1.0, r=0.03, z=4)
ax.plot([2.25, ZX0], [-0.02, ZY + 0.28], color=LIGHT, lw=0.5, ls='--'); ax.plot([2.95, ZX0 + ZW], [-0.02, ZY + 0.28], color=LIGHT, lw=0.5, ls='--')
for i in range(nres):
    ax.add_patch(Rectangle((ZX0 + i * ZW / nres, ZY), ZW / nres, 0.28, fc=PINK if dis[i] else GREEN, ec='white', lw=0.4, zorder=3))
ax.add_patch(Rectangle((ZX0, ZY), ZW, 0.28, fc='none', ec=INK, lw=0.6, zorder=4))
ax.plot([ZX0 + 6 * ZW / nres, ZX0 + 17 * ZW / nres], [ZY - 0.1] * 2, color=PINK_E, lw=1.2)
ax.text(ZX0 + 11.5 * ZW / nres, ZY - 0.22, 'IDR (≥20 disordered residues in the protein)', ha='center', va='top', fontsize=5.0, color=PINK_E)
ax.text(ZX0, ZY - 0.62, 'residue 1', fontsize=5.2, color=MUTED, va='center'); ax.text(ZX0 + ZW, ZY - 0.62, 'residue n', fontsize=5.2, color=MUTED, va='center', ha='right')
# definitions of the core metrics
DEF = [('IDR', 'run of ≥20 residues predicted disordered (MetaPredict v3)'),
       ('%IDR', 'percentage of a protein\'s residues inside IDRs'),
       ('exon IDR fraction', 'share of an exon\'s residues inside IDRs (0–1)'),
       ('constitutive / alt. exon', 'exon in every / only some of the gene\'s transcripts'),
       ('alt. − constitutive', 'difference in mean exon IDR fraction between\nalternative and constitutive exons'),
       ('shared / isoform-specific', 'codon encoded by every / only some isoforms'),
       ('IDR outcome', 'an IDR followed into another isoform: intact, truncated\n(ends lost), removed, or internally deleted'),
       ('gene age', 'time since the gene family arose (GenOrigin, Myr)'),
       ('splice-site conservation', 'mammals (of 405) keeping the exon\'s GT/AG;\npoorly conserved: fewer than 100')]
DX, DY = 5.0, 1.2
ax.text(DX, DY, 'Core metrics', fontsize=6.6, fontweight='bold', va='center')
y = DY - 0.3
for term, d in DEF:
    ax.text(DX, y, term, fontsize=5.3, fontweight='bold', va='top', color=INK)
    ax.text(DX + 1.8, y, d, fontsize=5.1, va='top', color=INK, linespacing=1.15)
    y -= 0.26 + 0.17 * d.count('\n')

# ---------------- b: approach
ax = fig.add_subplot(gs[1]); ax.set_xlim(0, 10); ax.set_ylim(0.05, 2.25); ax.axis('off'); tag(ax, 'b', 0.0, 0.95)
ax.text(0.35, 2.05, 'Approach', fontsize=8, fontweight='bold', va='center')
STEPS = [('Reviewed human\nisoforms', f'{CC["isoforms"]:,} isoforms\n{CC["genes"]:,} genes'),
         ('Predicted IDRs\n(MetaPredict)', f'{CC["idrs"]:,} IDRs\n≥20 residues'),
         ('Exact transcript\nmatch, GRCh38', f'{CC["main"]:,} isoforms\n{CC["analysis_unique_exons"]:,} exons'),
         ('Residues matched by\ngenomic codon', f'{RF["positions"] / 1e6:.1f} M positions\nall isoforms compared globally'),
         ('Classes, families, gene\nage, splice-site conservation', f'11 classes · {PC["I_n_families_ge15"]} families\n405 mammals')]
w = 1.72; xs = np.linspace(0.35, 10 - 0.35 - w, 5)
for k, ((t, sub), x) in enumerate(zip(STEPS, xs)):
    sc = ['#3f7fbf', DIS_E, '#2a9d8f', '#7b5ea7', '#c99a06'][k]
    ax.scatter([x + w / 2], [1.5], s=170, color=sc, zorder=3, lw=0); ax.text(x + w / 2, 1.5, str(k + 1), ha='center', va='center', color='white', fontsize=7, fontweight='bold', zorder=4)
    ax.text(x + w / 2, 0.95, t, ha='center', va='center', fontsize=6.2, fontweight='bold', color=sc)
    ax.text(x + w / 2, 0.38, sub, ha='center', va='center', fontsize=5.8, color=MUTED)
    if k < 4: arrow(ax, (x + w / 2 + 0.3, 1.5), (xs[k + 1] + w / 2 - 0.3, 1.5), c=LIGHT, ms=6)

# ---------------- c: findings
ax0 = fig.add_subplot(gs[2]); ax0.axis('off'); tag(ax0, 'c', 0.0, 1.2)
ax0.text(0.035, 1.23, 'Findings', fontsize=8, fontweight='bold', va='center', transform=ax0.transAxes)
sub = gs[2].subgridspec(1, 5, wspace=0.75)
FT = [('Disorder is set by\nclass and family', ), ('Alternative exons are more\ndisordered in every class', ), ('Splicing changes IDRs\nmostly at their ends', ),
      ('Disorder-to-order is a\nbetween-family pattern', ), ('Exons with few-mammal splice\nsites are more disordered', )]
for k in range(5):
    ax = fig.add_subplot(sub[0, k]); ax.set_title(FT[k][0], fontsize=6.1, fontweight='bold', pad=6)
    for s_ in ax.spines.values(): s_.set_visible(False)
    ax.tick_params(labelsize=5, length=2)
    if k == 0:
        v = [A.loc[c, 'canonical_median'] for c in CL]; ax.bar(range(len(v)), v, color=[plt.get_cmap('RdPu')(0.25 + 0.6 * i / len(v)) for i in range(len(v))], width=0.8)
        ax.set_xticks([0, len(v) - 1]); ax.set_xticklabels(['enzymes', 'DNA-binding\nTFs'], fontsize=4.8); ax.set_ylabel('median % disordered', fontsize=5.2)
        ax.text(0.03, 0.97, f'family explains {100 * PC["I_variance_explained_pct_idr"]["family"]:.0f}%\nof variance', transform=ax.transAxes, fontsize=4.8, va='top', color=MUTED)
    elif k == 1:
        sel = ['Metabolic enzymes', 'Receptors & channels', 'DNA-binding TFs']; x = np.arange(3)
        ax.bar(x - 0.2, [C.loc[c, 'pct_idr_exon_const'] for c in sel], 0.38, color=MUTED, label='constitutive')
        ax.bar(x + 0.2, [C.loc[c, 'pct_idr_exon_alt'] for c in sel], 0.38, color=ORANGE, label='alternative')
        ax.set_xticks(x); ax.set_xticklabels(['enzymes', 'receptors', 'TFs'], fontsize=4.6, rotation=30, ha='right'); ax.set_ylabel('% exons encoding an IDR', fontsize=5.2)
        ax.legend(frameon=False, fontsize=4.6, loc='upper left', handlelength=0.8)
        ax.text(0.03, 0.8, '3 of 11 classes shown;\nholds in all 11', transform=ax.transAxes, fontsize=4.6, va='top', color=MUTED)
    elif k == 2:
        P = RF['all_pairs_fate_pct_pooled']; FA = [('intact', LIGHT), ('truncated', AQUA), ('removed', ORANGE), ('internal deletion', BLUE)]
        b = 0
        for f, c in FA:
            ax.bar(0, P[f], bottom=b, color=c, width=0.6, edgecolor='white', lw=0.5); ax.text(0.38, b + P[f] / 2, f'{f} {P[f]:.0f}%', fontsize=4.8, va='center'); b += P[f]
        ax.set_xlim(-0.4, 2.0); ax.set_ylim(0, 100); ax.set_xticks([]); ax.set_ylabel('IDRs followed between\nisoforms (%)', fontsize=5.2)
    elif k == 3:
        S = EV['E6_slope_pct_idr_per_100Myr']; v = [S['all genes'], S['family genes, family fixed effects']]
        ax.bar([0, 1], [v[0][0], v[1][0]], color=[BLUE, ORANGE], width=0.6)
        ax.errorbar([0, 1], [v[0][0], v[1][0]], yerr=[[v[0][0] - v[0][1], v[1][0] - v[1][1]], [v[0][2] - v[0][0], v[1][2] - v[1][0]]], fmt='none', ecolor=INK, lw=0.6)
        refline(ax, 0, 'y'); ax.set_xticks([0, 1]); ax.set_xticklabels(['all genes', 'within\nfamilies'], fontsize=4.8)
        ax.set_ylabel('Δ %IDR per 100 Myr\nolder', fontsize=5.2); ax.set_ylim(-0.7, 0.25)
    else:
        X1 = EA['X1']; lab = list(X1)
        ax.bar([0, 1], [X1[lab[0]]['mean_idr'], X1[lab[2]]['mean_idr']], color=[ORANGE, MUTED], width=0.6)
        ax.set_xticks([0, 1]); ax.set_xticklabels(['<100\ngenomes', '≥350\ngenomes'], fontsize=4.8); ax.set_ylabel('mean exon IDR fraction', fontsize=5.2); ax.set_ylim(0, 0.56)
        ax.text(0.5, 0.97, f'within genes +{EA["X2_fe_young"]["coef"]:.2f}\nacidic −{abs(EA["X4_young_minus_old_disorder_composition"]["acid"]["diff_pp"]):.1f} pp',
                transform=ax.transAxes, fontsize=4.8, va='top', ha='center', color=MUTED)
save(fig, 'main/fig1_graphical_abstract')
