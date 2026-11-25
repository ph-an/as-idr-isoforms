"""Main Figure 6: disorder and splicing through gene age, by protein class.

a, %IDR of every isoform (canonical and alternative) by gene age (GenOrigin), for all genes and each protein class; violins
over isoforms, white line at the median; rho, Spearman correlation of gene age with the gene's mean %IDR over its isoforms
(one value per gene); b, length of the longest IDR of every isoform by gene age (median, line; interquartile range, band;
isoforms without an IDR count as 0); c, share of each gene's unique exons that are alternative, by gene age; d, alternative
minus constitutive exon IDR fraction by gene age (95% intervals from resampling genes); e, IDR outcomes across all ordered
isoform pairs by gene age. Family-level versions of a and b are Extended Data Fig. 8; splice-site conservation is Fig. 7.
"""
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from style import *

EVd = AN / 'evolution/outputs'
PC = json.load(open(AN / 'protein_classes/outputs/key_numbers.json')); CL = PC['classes_by_median_pct_idr']
X = pd.read_csv(EVd / 'isoform_disorder.csv'); B = pd.read_csv(EVd / 'isoform_disorder_by_age.csv'); T = pd.read_csv(EVd / 'isoform_disorder_trend.csv')
E2 = pd.read_csv(EVd / 'E2_splicing_by_age.csv').set_index('age_bin').reindex(AGE)
AGEL = ['<100', '100–500', '500–1,000', '>1,000']; SH = ['#c6dbef', '#6baed6', '#2171b5', '#08306b']
GROUPS = ['All'] + CL[::-1]
NAME = {'All': 'all genes', **{c: SHORT[c] for c in CL}, 'Transcription co-regulators': 'TF co-regulators', 'Chromatin & RNA processing': 'Chromatin & RNA',
        'Cytoskeletal & structural': 'Cytoskeletal', 'Other / unclassified': 'Other / unclassified', 'Kinases & phosphatases': 'Kinases & phosph.'}


def rho_text(ax, g, m, y=-0.1, va='top'):
    t = T[(T.group == g) & (T.kind != 'family')].iloc[0]; r, p = t[f'{m}_rho'], t[f'{m}_p']
    c = (ORANGE if r > 0 else BLUE) if p < 0.05 else MUTED
    ax.text(0.5, y, f'ρ = {r:+.2f}', transform=ax.transAxes, ha='center', va=va, fontsize=5.4, color=c)


fig = plt.figure(figsize=(7.1, 9.3)); G = fig.add_gridspec(3, 1, height_ratios=[1.2, 1.15, 1.0], hspace=0.36, top=0.93, bottom=0.05)
# a: %IDR of all isoforms, violins per class
ga = G[0].subgridspec(2, 6, hspace=0.5, wspace=0.18)
for k, g in enumerate(GROUPS):
    ax = fig.add_subplot(ga[k // 6, k % 6]); z = X if g == 'All' else X[X.pclass == g]
    for j, b in enumerate(AGE):
        v = z[z.age_bin == b].pct_idr.values
        if len(v) >= 20:
            vp = ax.violinplot(v, [j], widths=0.85, showextrema=False)
            for body in vp['bodies']: body.set_facecolor(SH[j]); body.set_edgecolor('none'); body.set_alpha(1)
            ax.plot([j - 0.25, j + 0.25], [np.median(v)] * 2, color='white' if j >= 2 else INK, lw=1.2)
        elif len(v): ax.plot(j, np.median(v), 'o', ms=2.5, color=SH[j])
    ax.set_ylim(0, 100); ax.set_xlim(-0.6, 3.6); ax.set_xticks([]); grid(ax)
    ax.set_yticks([0, 50, 100]); ax.tick_params(axis='y', labelsize=5.2)
    if k % 6: ax.set_yticklabels([])
    else: ax.set_ylabel('%IDR (all isoforms)', fontsize=5.8)
    ax.set_title(NAME[g], fontsize=5.6, pad=2, fontweight='bold' if g == 'All' else 'normal'); rho_text(ax, g, 'pct_idr')
    if k == 0: tag(ax, 'a', -0.55, 1.2)
fig.legend([Patch(color=c) for c in SH], [a + ' Myr' for a in AGEL], ncol=4, frameon=False, fontsize=5.8, loc='upper center', bbox_to_anchor=(0.5, 0.985), handlelength=1, title='gene age (GenOrigin); ρ, Spearman correlation over genes (blue, negative; orange, positive; grey, P ≥ 0.05)', title_fontsize=5.8)
# b: longest IDR per isoform, median and IQR
gb = G[1].subgridspec(2, 6, hspace=0.62, wspace=0.18)
for k, g in enumerate(GROUPS):
    ax = fig.add_subplot(gb[k // 6, k % 6]); z = B[(B.group == g) & (B.kind != 'family')].set_index('age_bin').reindex(AGE)
    ok = z.isoforms.values >= 20; x = np.arange(4)
    ax.fill_between(x[ok], z.longest_q1[ok], z.longest_q3[ok], color=BLUE, alpha=0.15, lw=0)
    ax.plot(x[ok], z.longest_median[ok], '-', color=BLUE, lw=1.1)
    ax.scatter(x[ok], z.longest_median[ok], s=9, c=[SH[i] for i in x[ok]], edgecolors=BLUE, linewidths=0.5, zorder=3)
    ax.set_ylim(0, 450); ax.set_xlim(-0.4, 3.4); grid(ax); ax.set_yticks([0, 200, 400]); ax.tick_params(axis='y', labelsize=5.2)
    ax.set_xticks(x); ax.set_xticklabels(['<100', '', '', '>1,000'] if k // 6 else [], fontsize=4.8)
    if k == 6: ax.text(0, -0.42, 'gene age (Myr) →', transform=ax.transAxes, fontsize=5.2, color=MUTED)
    if k % 6: ax.set_yticklabels([])
    else: ax.set_ylabel('longest IDR (residues)', fontsize=5.8)
    ax.set_title(NAME[g], fontsize=5.6, pad=2, fontweight='bold' if g == 'All' else 'normal'); rho_text(ax, g, 'longest', y=0.97)
    if k == 0: tag(ax, 'b', -0.55, 1.2)
# c, d, e
gc = G[2].subgridspec(1, 3, wspace=0.75, width_ratios=[1, 1, 1.15]); x = np.arange(4)
ax = fig.add_subplot(gc[0]); ax.plot(x, 100 * E2.alt_exon_fraction, 'o-', color=INK, lw=1.3, ms=4)
for xi, v in zip(x, E2.alt_exon_fraction): ax.text(xi, 100 * v + 2, f'{100 * v:.0f}%', ha='center', fontsize=5.4)
ax.set_ylim(40, 80); ax.set_ylabel("alternative share of a gene's\nunique exons (%)"); ax.set_xticks(x); ax.set_xticklabels(AGEL, rotation=30, ha='right', fontsize=5.8)
ax.set_xlabel('gene age (Myr)'); grid(ax); tag(ax, 'c', -0.45)
ax = fig.add_subplot(gc[1])
ax.errorbar(x, E2.coupling, yerr=[E2.coupling - E2.coupling_lo, E2.coupling_hi - E2.coupling], fmt='s', color=ORANGE, ms=4, lw=1, capsize=0)
ax.set_ylim(0, 0.25); ax.set_ylabel('alternative − constitutive\nexon IDR fraction'); ax.set_xticks(x); ax.set_xticklabels(AGEL, rotation=30, ha='right', fontsize=5.8)
ax.set_xlabel('gene age (Myr)'); refline(ax, 0, 'y', 'no difference', ha='left', x=0.02); grid(ax); tag(ax, 'd', -0.45)
ax = fig.add_subplot(gc[2]); FATES = [('intact', LIGHT), ('internal deletion', BLUE), ('truncated', AQUA), ('removed', ORANGE)]
left = np.zeros(4); yb = x[::-1]
for f, c in FATES:
    v = E2[f'fate_{f}'].values; ax.barh(yb, v, left=left, color=c, height=0.7, edgecolor='white', lw=0.5)
    for yi, l, vi in zip(yb, left, v):
        if vi >= 2.5: ax.text(l + vi / 2, yi, f'{vi:.0f}', ha='center', va='center', fontsize=5.4 if vi >= 8 else (4.6 if vi >= 5 else 3.8), color=INK if f in ('intact', 'truncated') else 'white')
    left += v
ax.set_yticks(yb); ax.set_yticklabels([a + ' Myr' for a in AGEL], fontsize=5.6); ax.set_xlim(0, 100); ax.set_xlabel('IDR outcomes across\nisoform pairs (%)')
ax.legend([Patch(color=c) for _, c in FATES], [f for f, _ in FATES], ncol=2, frameon=False, fontsize=5.2, loc='lower left', bbox_to_anchor=(-0.35, 1.0), handlelength=1, columnspacing=0.6)
tag(ax, 'e', -0.42, 1.25)
save(fig, 'main/fig6_evolution')
