"""Extended Data Figure 8: disorder of all isoforms through gene age within protein families.

The 15 UniProt families with at least 15 dated genes and at least five genes in each of three or more age bins. a, %IDR of
every isoform by gene age (violins over isoforms; white line, median); b, length of the longest IDR of every isoform (median,
line; interquartile range, band). Age bins with fewer than five genes are left empty. rho, Spearman correlation of gene age
with the gene's mean over its isoforms (blue, negative; orange, positive; grey, P >= 0.05). Numbers from
analysis/evolution/isoform_disorder_by_age.py.
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from style import *

EVd = AN / 'evolution/outputs'
X = pd.read_csv(EVd / 'isoform_disorder.csv'); B = pd.read_csv(EVd / 'isoform_disorder_by_age.csv'); T = pd.read_csv(EVd / 'isoform_disorder_trend.csv')
F = B[B.kind == 'family']; ok = F[F.genes >= 5].groupby('group').size(); fams = ok[ok >= 3].index
T = T[(T.kind == 'family') & T.group.isin(fams)].sort_values('pct_idr_rho', ascending=False); FAMS = list(T.group)
SH = ['#c6dbef', '#6baed6', '#2171b5', '#08306b']; AGEL = ['<100', '100–500', '500–1,000', '>1,000']
short = lambda f: f.replace(' superfamily', ' (super)').replace(' family', '').replace('Krueppel C2H2-type zinc-finger protein', 'C2H2 zinc finger').replace('TRAFAC class myosin-kinesin ATPase', 'Myosin-kinesin ATPase').replace(' (TC 2.A.29)', '').replace('G-protein coupled receptor 1', 'GPCR 1')


def rho(ax, f, m, y):
    t = T[T.group == f].iloc[0]; r, p = t[f'{m}_rho'], t[f'{m}_p']; c = (ORANGE if r > 0 else BLUE) if p < 0.05 else MUTED
    ax.text(0.5, y, f'ρ = {r:+.2f}', transform=ax.transAxes, ha='center', va='top', fontsize=5.2, color=c)


fig = plt.figure(figsize=(7.1, 9.2)); G = fig.add_gridspec(2, 1, hspace=0.22, top=0.94, bottom=0.04)
ga = G[0].subgridspec(3, 5, hspace=0.6, wspace=0.18); gb = G[1].subgridspec(3, 5, hspace=0.6, wspace=0.18)
for k, f in enumerate(FAMS):
    z = X[X.family == f]; ng = z.groupby('age_bin').base_accession.nunique()
    ax = fig.add_subplot(ga[k // 5, k % 5])
    for j, b in enumerate(AGE):
        if ng.get(b, 0) < 5: continue
        v = z[z.age_bin == b].pct_idr.values
        vp = ax.violinplot(v, [j], widths=0.85, showextrema=False)
        for body in vp['bodies']: body.set_facecolor(SH[j]); body.set_edgecolor('none'); body.set_alpha(1)
        ax.plot([j - 0.25, j + 0.25], [np.median(v)] * 2, color='white' if j >= 2 else INK, lw=1.1)
    ax.set_ylim(0, 100); ax.set_xlim(-0.6, 3.6); ax.set_xticks([]); grid(ax); ax.set_yticks([0, 50, 100]); ax.tick_params(axis='y', labelsize=5)
    if k % 5: ax.set_yticklabels([])
    else: ax.set_ylabel('%IDR (all isoforms)', fontsize=5.6)
    ax.set_title(f'{short(f)} ({int(T[T.group == f].genes.iloc[0])})', fontsize=5.3, pad=2); rho(ax, f, 'pct_idr', -0.06)
    if k == 0: tag(ax, 'a', -0.55, 1.2)
    ax = fig.add_subplot(gb[k // 5, k % 5]); zz = F[F.group == f].set_index('age_bin').reindex(AGE); m = (zz.genes >= 5).values; x = np.arange(4)
    ax.fill_between(x[m], zz.longest_q1[m], zz.longest_q3[m], color=BLUE, alpha=0.15, lw=0); ax.plot(x[m], zz.longest_median[m], '-', color=BLUE, lw=1)
    ax.scatter(x[m], zz.longest_median[m], s=8, c=[SH[i] for i in x[m]], edgecolors=BLUE, linewidths=0.5, zorder=3)
    ax.set_ylim(0, 450); ax.set_xlim(-0.4, 3.4); grid(ax); ax.set_yticks([0, 200, 400]); ax.tick_params(axis='y', labelsize=5)
    ax.set_xticks(x); ax.set_xticklabels(['<100', '', '', '>1,000'] if k // 5 == 2 else [], fontsize=4.8)
    if k % 5: ax.set_yticklabels([])
    else: ax.set_ylabel('longest IDR (residues)', fontsize=5.6)
    ax.set_title(short(f), fontsize=5.3, pad=2); rho(ax, f, 'longest', 0.97)
    if k == 0: tag(ax, 'b', -0.55, 1.2)
fig.legend([Patch(color=c) for c in SH], [a + ' Myr' for a in AGEL], ncol=4, frameon=False, fontsize=5.8, loc='upper center', bbox_to_anchor=(0.5, 0.995), handlelength=1, title='gene age (GenOrigin); families ordered by ρ for %IDR; genes in parentheses', title_fontsize=5.8)
save(fig, 'extended_data/ed8_family_isoform_age')
