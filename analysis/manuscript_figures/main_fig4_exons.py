"""Main Figure 4: alternative exons encode more disorder than constitutive exons in every protein class; IDRs are modestly exon-contained.

a, % of unique exons encoding an IDR (IDR fraction >= 0.5), alternative vs constitutive, per class; b, per-gene
difference in mean exon IDR fraction (alternative - constitutive); c, the class-level difference with gene-bootstrap 95%
intervals, and the all-gene estimate under six estimators; d, share of each gene's unique exons that are alternative;
e, exon boundaries inside IDR and ordered segments relative to a
length-matched placement null, by segment length (lines) and by protein class (points).
"""
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from style import *

PC = AN / 'protein_classes' / 'outputs'
K = json.load(open(PC / 'key_numbers.json')); CL = K['classes_by_median_pct_idr']
TC = pd.read_csv(PC / 'C_exon_coupling_by_class.csv', index_col=0).loc[CL]
D = pd.read_csv(PC / 'D_gene_delta.csv', index_col=0); E = pd.read_csv(PC / 'E_gene_alt_fraction.csv', index_col=0)
F3 = json.load(open(AN / 'manuscript_figures/outputs/fig3_numbers.json'))
RF = json.load(open(AN / 'reference_free/outputs/key_numbers.json'))
fig = plt.figure(figsize=(7.1, 6.3)); gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.45, width_ratios=[1.25, 1, 1])
y = np.arange(len(CL))
lab = [f'{SHORT[c]} ({int(TC.loc[c, "genes"]):,})' for c in CL]
# a
ax = fig.add_subplot(gs[0, 0])
ax.barh(y + 0.2, TC.pct_idr_exon_const, 0.38, color=MUTED, label='constitutive')
ax.barh(y - 0.2, TC.pct_idr_exon_alt, 0.38, color=ORANGE, label='alternative')
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=6); ax.set_xlabel('exons encoding an IDR (%)'); grid(ax, 'x')
ax.legend(frameon=False, fontsize=6, loc='lower right'); tag(ax, 'a', -0.75)
# b
ax = fig.add_subplot(gs[0, 1])
for yi, c in zip(y, CL):
    v = D[D.pclass == c].delta.values; q = np.percentile(v, [10, 25, 50, 75, 90])
    ax.plot([q[0], q[4]], [yi, yi], color=MUTED, lw=0.7); ax.add_patch(plt.Rectangle((q[1], yi - 0.3), q[3] - q[1], 0.6, color=BLUE, alpha=0.7, lw=0))
    ax.plot([q[2]] * 2, [yi - 0.3, yi + 0.3], color=INK, lw=1)
    ax.text(1.02, yi, f'{100 * TC.loc[c, "frac_genes_alt_higher"]:.0f}%', va='center', fontsize=5.6, color=MUTED, transform=ax.get_yaxis_transform())
refline(ax, 0, 'x', 'no difference'); ax.set_yticks(y); ax.set_yticklabels([]); ax.set_xlim(-1, 1)
ax.set_xlabel('per gene: alternative − constitutive\nmean exon IDR fraction'); grid(ax, 'x'); tag(ax, 'b', -0.05)
ax.text(1.02, len(CL) - 0.2, 'genes\n>0', fontsize=5.6, color=MUTED, transform=ax.get_yaxis_transform(), va='bottom')
# c
ax = fig.add_subplot(gs[0, 2])
ax.hlines(y, TC.diff_lo, TC.diff_hi, color=INK, lw=1); ax.plot(TC['diff'], y, 'o', ms=4, color=ORANGE, mec='white', mew=0.5)
ax.axvline(F3['boot'][0], color=ORANGE, lw=0.6, ls='--'); refline(ax, 0, 'x', 'no\ndifference', ha='left')
ax.set_yticks(y); ax.set_yticklabels([]); ax.set_xlim(-0.02, 0.3); ax.set_xlabel('alternative − constitutive\nexon IDR fraction (95% CI)'); grid(ax, 'x'); tag(ax, 'c', -0.05)
ax.text(F3['boot'][0], len(CL) - 0.3, ' all genes', fontsize=5.6, color=ORANGE, va='bottom')
# d
ax = fig.add_subplot(gs[1, 0])
for yi, c in zip(y, CL):
    v = E[E.pclass == c].alt_fraction.values; q = np.percentile(v, [10, 25, 50, 75, 90])
    ax.plot([q[0], q[4]], [yi, yi], color=MUTED, lw=0.7); ax.add_patch(plt.Rectangle((q[1], yi - 0.3), q[3] - q[1], 0.6, color=AQUA, alpha=0.7, lw=0))
    ax.plot([q[2]] * 2, [yi - 0.3, yi + 0.3], color=INK, lw=1)
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=6); ax.set_xlabel("share of a gene's unique exons\nthat are alternative"); ax.set_xlim(0, 1); grid(ax, 'x'); tag(ax, 'd', -0.75)
# f: boundaries
R6 = pd.read_csv(AN / 'reanalysis_2026-10-26/outputs/R6_by_length.csv'); FO = pd.read_csv(PC / 'F_modularity_by_class.csv')
r1 = gs[1, 1:].subgridspec(1, 2, wspace=0.12, width_ratios=[1.1, 0.9])
ax = fig.add_subplot(r1[0]); bins = ['20-50', '51-100', '101-200', '201-400', '>400']
for kind, c in [('IDR', BLUE), ('ordered', MUTED)]:
    z = R6[R6.kind == kind].set_index('lbin').reindex(bins); ax.plot(range(5), z.O_E, 'o-', color=c, lw=1.5, ms=4, label=f'{kind} segments')
refline(ax, 1, 'y', 'random expectation', ha='left', x=0.02); ax.set_xticks(range(5)); ax.set_xticklabels([b.replace('-', '–') for b in bins], rotation=35, ha='right', fontsize=5.8)
ax.set_xlabel('segment length (residues)'); ax.set_ylabel('exon boundaries in segment\n(observed / expected)', fontsize=6); ax.legend(frameon=False, fontsize=6); grid(ax); tag(ax, 'e', -0.38)
ax = fig.add_subplot(r1[1]); FI = FO[FO.kind == 'IDR'].set_index('pclass').reindex(CL); FOr = FO[FO.kind == 'ordered'].set_index('pclass').reindex(CL)
ax.plot(FOr.O_E, y, 's', ms=3, color=MUTED, label='ordered'); ax.plot(FI.O_E, y, 'o', ms=3.5, color=BLUE, label='IDR')
refline(ax, 1, 'x', 'random\nexpectation', ha='left'); ax.set_yticks(y); ax.set_yticklabels([SHORT[c] for c in CL], fontsize=5.6); ax.yaxis.tick_right()
ax.set_xlabel('observed / expected\nboundaries, by class'); ax.legend(frameon=False, fontsize=5.6, loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2, handlelength=1); grid(ax, 'x')
save(fig, 'main/fig4_exons')
