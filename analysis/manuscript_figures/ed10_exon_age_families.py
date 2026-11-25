"""Extended Data Figure 10: exon IDR fraction by splice-site conservation within UniProt families.

Families with >= 15 genes and >= 10 internal exons in every conservation bin, ordered by the within-family difference.
a, mean exon IDR fraction by the number of 405 mammalian genomes keeping the exon's less conserved splice site (number of
exons in small type); b, mean IDR fraction of exons conserved in < 350 genomes minus exons conserved in >= 380 genomes,
gene-bootstrap 95% CI; filled, interval excludes 0; colour, the family's most common protein class.
"""
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from style import *

X = pd.read_csv(AN / 'exon_age/outputs/X6_exon_idr_by_conservation_by_family.csv').sort_values('diff').reset_index(drop=True)
K = json.load(open(AN / 'exon_age/outputs/key_numbers.json'))
B = ['<350', '350-379', '380-395', '396-405']
CLS = sorted(X.pclass.unique()); pal = plt.get_cmap('tab20'); col = {c: pal(i / max(1, len(CLS) - 1)) for i, c in enumerate(CLS)}
fig = plt.figure(figsize=(7.1, 6.8)); gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.8], wspace=0.06)
ax = fig.add_subplot(gs[0, 0]); y = np.arange(len(X))
M = X[[f'mean_{b}' for b in B]].values; N = X[[f'n_{b}' for b in B]].values
im = ax.imshow(M, cmap='Blues', vmin=0, vmax=0.8, aspect='auto', origin='lower')
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i + 0.08, f'{M[i, j]:.2f}', ha='center', va='center', fontsize=5.2, color='white' if M[i, j] > 0.45 else INK)
        ax.text(j + 0.45, i - 0.32, f'{N[i, j]}', ha='right', va='bottom', fontsize=3.6, color='white' if M[i, j] > 0.45 else MUTED)
lab = [f'{f.replace(" family", "").replace(" superfamily", " (super)").replace("Krueppel C2H2-type zinc-finger protein", "C2H2 zinc finger")[:42]} ({g})' for f, g in zip(X.family, X.genes)]
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=5.4); ax.set_xticks(range(4)); ax.set_xticklabels(['<350', '350–379', '380–395', '396–405'], fontsize=6)
ax.set_xlabel('mammals (of 405) keeping the exon’s\nless conserved splice site'); [s.set_visible(False) for s in ax.spines.values()]
cax = ax.inset_axes([0.15, 1.035, 0.7, 0.018]); cb = fig.colorbar(im, cax=cax, orientation='horizontal'); cax.xaxis.set_ticks_position('top'); cax.xaxis.set_label_position('top'); cb.set_label('mean exon IDR fraction (internal exons)', fontsize=6); cb.ax.tick_params(labelsize=5.5)
tag(ax, 'a', -0.75, 1.03)
ax = fig.add_subplot(gs[0, 1], sharey=ax)
for yi, r in zip(y, X.itertuples()):
    sig = r.lo > 0 or r.hi < 0
    ax.plot([r.lo, r.hi], [yi, yi], color=MUTED, lw=0.7); ax.plot(r.diff, yi, 'o', ms=3.8, color=col[r.pclass] if sig else 'white', mec=col[r.pclass], mew=1)
refline(ax, 0, 'x', 'no difference', ha='right'); ax.axvline(K['X6_median_diff'], color=ORANGE, lw=0.6, ls='--')
ax.tick_params(axis='y', labelleft=False); ax.set_xlabel('less (<350) minus more (≥380)\nconserved exons, IDR fraction (95% CI)'); grid(ax, 'x')
ax.text(K['X6_median_diff'], len(X) - 0.2, ' median', color=ORANGE, fontsize=5.4, va='bottom')
ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', color=col[c], ms=3.8, label=SHORT.get(c, c)) for c in CLS], frameon=False, fontsize=5.2, loc='upper center', bbox_to_anchor=(0.5, -0.11), ncol=3, title='most common class', title_fontsize=5.4, columnspacing=0.8, handletextpad=0.3)
tag(ax, 'b', -0.05, 1.03)
save(fig, 'extended_data/ed10_exon_age_families')
