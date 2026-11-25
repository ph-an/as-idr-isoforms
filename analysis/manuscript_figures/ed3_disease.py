"""Extended Data Figure 3: disease inheritance and disorder.

Genes with autosomal dominant, autosomal recessive or both modes of disease inheritance in UniProt ('Involvement in
disease'), and genes with no disease annotation: a, median IDR fraction of the canonical isoform; b, alternative minus
constitutive exon IDR fraction; c, isoform-specific minus shared IDR fraction. Bars, 95% intervals from resampling genes.
Numbers from analysis/evolution/outputs/E4_inheritance.csv (formerly Fig. 3c).
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from style import *

E4 = pd.read_csv(AN / 'evolution/outputs/E4_inheritance.csv').set_index('group')
grp = ['AD only', 'AR only', 'AD and AR', 'no disease text']; cols = [ORANGE, BLUE, MUTED, LIGHT]
names = {'AD only': 'dominant only', 'AR only': 'recessive only', 'AD and AR': 'dominant and recessive', 'no disease text': 'no disease annotation'}
MS = [('median_pct_idr', 'median_lo', 'median_hi', 100, 'median IDR fraction\n(canonical isoform)'),
      ('coupling', 'lo', 'hi', 1, 'alternative − constitutive\nexon IDR fraction'),
      ('var_minus_core', 'vc_lo', 'vc_hi', 1, 'isoform-specific −\nshared IDR fraction')]
fig = plt.figure(figsize=(7.1, 2.7)); gs = fig.add_gridspec(1, 3, wspace=0.55, top=0.8, bottom=0.18)
for k, (m, lo, hi, sc, yl) in enumerate(MS):
    ax = fig.add_subplot(gs[0, k]); x = np.arange(4)
    v = [E4.loc[g, m] / sc for g in grp]; l = [E4.loc[g, lo] / sc for g in grp]; h = [E4.loc[g, hi] / sc for g in grp]
    ax.bar(x, v, 0.75, color=cols); ax.errorbar(x, v, yerr=[np.subtract(v, l), np.subtract(h, v)], fmt='none', ecolor=INK, lw=0.7, capsize=1.5)
    ax.set_xticks(x); ax.set_xticklabels(['dominant', 'recessive', 'both', 'none'], fontsize=5.8, rotation=30, ha='right'); ax.set_ylabel(yl)
    ax.set_ylim(0, 0.4 if k == 0 else 0.27); grid(ax); tag(ax, 'abc'[k], -0.32)
fig.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in cols], [f'{names[g]} ({int(E4.loc[g, "genes"]):,} genes)' for g in grp], frameon=False, fontsize=6,
           loc='upper center', bbox_to_anchor=(0.5, 0.99), ncol=4, handlelength=1)
save(fig, 'extended_data/ed3_disease')
