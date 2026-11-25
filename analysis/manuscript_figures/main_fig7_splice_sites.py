"""Main Figure 7: exons with poorly conserved splice sites encode more, and different, disorder.

a, mean IDR fraction of internal exons (bars, 95% intervals; left axis) and share of them that are alternative (line; right
axis) by the number of 405 mammalian genomes that keep the GT/AG of the exon's less conserved splice site; b, composition of disorder from poorly conserved minus deeply conserved exons of the same
genes, beside the shared disorder of young minus old genes.
"""
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from style import *

EV = json.load(open(AN / 'evolution/outputs/key_numbers.json')); EA = json.load(open(AN / 'exon_age/outputs/key_numbers.json'))
D = pd.read_csv(AN / 'exon_age/outputs/exon_conservation.csv'); I = D[(D.pos == 'internal') & (D.n_sites == 2)]
R = pd.read_csv(AN / 'exon_age/outputs/robustness_models.csv').set_index('model'); ER = json.load(open(AN / 'exon_age/outputs/robustness.json'))
S = json.load(open(AN / 'sensitivity/outputs/sensitivity.json'))
fig = plt.figure(figsize=(7.1, 3.1)); Gt = fig.add_gridspec(1, 1, left=0.07, right=0.55, top=0.8, bottom=0.25)
Gb = fig.add_gridspec(1, 1, left=0.74, right=0.99, top=0.8, bottom=0.25)
# a
ax = fig.add_subplot(Gt[0, 0])
edges = [-1, 49, 99, 199, 299, 349, 379, 395, 405]; lab = ['<50', '50–99', '100–199', '200–299', '300–349', '350–379', '380–395', '396–405']
g = I.assign(cd=pd.cut(I.cons, edges, labels=lab)).groupby('cd', observed=True).agg(m=('idr_frac', 'mean'), se=('idr_frac', lambda v: v.std() / np.sqrt(len(v))), alt=('alt', 'mean'), n=('idr_frac', 'size'))
xx = np.arange(len(g)); ax.bar(xx, g.m, yerr=1.96 * g.se, color=BLUE, width=0.62, error_kw=dict(lw=0.6, ecolor=INK)); ax.set_ylim(0, 0.5)
for xi, n in zip(xx, g.n): ax.text(xi, 0.012, f'{n:,}', ha='center', va='bottom', fontsize=4.8, color='white', rotation=90)
ax.set_xticks(xx); ax.set_xticklabels(lab, fontsize=5.4, rotation=35, ha='right'); ax.set_ylabel('mean exon IDR fraction\n(bars)', color=BLUE)
ax.set_xlabel('mammalian genomes (of 405) keeping the\nexon’s less conserved splice site')
ax.text(0.0, -0.36, '← poorly conserved', transform=ax.transAxes, fontsize=5.4, color=MUTED); ax.text(1.0, -0.36, 'deeply conserved →', transform=ax.transAxes, fontsize=5.4, color=MUTED, ha='right'); grid(ax); tag(ax, 'a', -0.12, 1.08)
ax2 = ax.twinx(); ax2.spines['right'].set_visible(True); ax2.plot(xx, 100 * g.alt, 'o-', color=ORANGE, ms=3.5, lw=1.3)
ax2.set_ylim(0, 100); ax2.set_ylabel('alternative exons\n(%, line)', color=ORANGE); ax2.tick_params(axis='y', colors=ORANGE)
# c
ax = fig.add_subplot(Gb[0, 0]); X4 = EA['X4_young_minus_old_disorder_composition']; EVc = EV['E3_young_minus_old_core_pp']
NM = {'acid': 'acidic', 'basic': 'basic', 'his': 'His', 'cys': 'Cys', 'aro': 'aromatic', 'hyd': 'aliphatic', 'pro': 'Pro', 'st': 'Ser, Thr', 'qn': 'Gln, Asn', 'gly': 'Gly', 'ncpr': 'net charge'}
ks = sorted(NM, key=lambda k: X4[k]['diff_pp']); y = np.arange(len(ks))
ax.barh(y + 0.18, [X4[k]['diff_pp'] for k in ks], 0.36, color=[ORANGE if X4[k]['q'] < 0.05 else LIGHT for k in ks])
ax.errorbar([X4[k]['diff_pp'] for k in ks], y + 0.18, xerr=[[X4[k]['diff_pp'] - X4[k]['lo'] for k in ks], [X4[k]['hi'] - X4[k]['diff_pp'] for k in ks]], fmt='none', ecolor=INK, lw=0.5)
ax.barh(y - 0.18, [EVc[k] for k in ks], 0.36, color=BLUE, alpha=0.6)
refline(ax, 0, 'x', ' no difference', ha='left', y=0.005); ax.set_yticks(y); ax.set_yticklabels([NM[k] for k in ks], fontsize=5.6)
ax.set_xlabel('young − old disorder\n(percentage points)'); grid(ax, 'x')
s = EA['X4_similarity_spearman']['vs_young_gene_core']
ax.legend([Patch(color=ORANGE), Patch(color=LIGHT), Patch(color=BLUE, alpha=0.6)], ['poorly − deeply conserved exons, same gene (q < 0.05)', 'same, not significant', 'young − old genes (<500 vs >1,000 Myr)'],
          frameon=False, fontsize=5.0, loc='lower left', bbox_to_anchor=(-0.55, 1.0), handlelength=1)
ax.text(0.03, 0.97, f'profiles ρ = {s[0]:.2f}', transform=ax.transAxes, ha='left', va='top', fontsize=5.4); tag(ax, 'b', -0.45, 1.3)
save(fig, 'main/fig7_splice_sites')
