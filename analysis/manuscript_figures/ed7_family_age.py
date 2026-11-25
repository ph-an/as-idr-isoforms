"""Extended Data Figure 7: gene age and disorder within UniProt families.

Spearman correlation of canonical %IDR with gene age within each UniProt family that has at least 15 dated genes in at least
two age bins (95% CI, Fisher z); filled, Benjamini-Hochberg q < 0.05; colour, the family's most common protein class.
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from style import *

T = pd.read_csv(AN / 'evolution/outputs/E6_age_disorder_by_family.csv').sort_values('rho').reset_index(drop=True)
CLS = sorted(T.pclass.unique()); pal = plt.get_cmap('tab20')
col = {c: pal(i / max(1, len(CLS) - 1)) for i, c in enumerate(CLS)}
fig, ax = plt.subplots(figsize=(7.1, 7.6))
y = np.arange(len(T))
for yi, r in zip(y, T.itertuples()):
    z = np.arctanh(np.clip(r.rho, -0.999, 0.999)); se = 1 / np.sqrt(max(r.genes - 3, 1)); lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    ax.plot([lo, hi], [yi, yi], color=MUTED, lw=0.7)
    ax.plot(r.rho, yi, 'o', ms=4, color=col[r.pclass] if r.q < 0.05 else 'white', mec=col[r.pclass], mew=1)
lab = [f'{f.replace(" family", "").replace(" superfamily", " (super)").replace("Krueppel C2H2-type zinc-finger protein", "C2H2 zinc finger")[:44]} ({g})' for f, g in zip(T.family, T.genes)]
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=5.4); refline(ax, 0, 'x', 'no correlation with gene age', fs=6); ax.set_ylim(-1, len(T))
ax.set_xlabel('Spearman ρ, canonical %IDR vs gene age within family (95% CI)'); grid(ax, 'x')
ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', color=col[c], ms=4, label=SHORT.get(c, c)) for c in CLS], frameon=False, fontsize=5.6, loc='lower right', title='most common class', title_fontsize=5.8)
save(fig, 'extended_data/ed7_family_age')
