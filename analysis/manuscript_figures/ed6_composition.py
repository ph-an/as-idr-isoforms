"""Extended Data Figure 6: the sequence of variable disorder depends on how often and where it is included.

Composition of disordered positions that are present in only some of a gene's isoforms, compared with the gene's core
disorder (positions in every isoform), split by inclusion (most, half or a minority of isoforms) and by position
(internal: flanked by core sequence on both sides; terminal: at an isoform end). a, internal variable disorder;
b, terminal variable disorder. Colour: difference in percentage points (variable - core, gene-paired); dot: BH q < 0.05.
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from style import *

C = pd.read_csv(AN / 'reference_free/outputs/composition_vs_core.csv')
FEAT = [('acid', 'acidic (D, E)'), ('ncpr', 'net charge per residue'), ('basic', 'basic (K, R)'), ('his', 'histidine'), ('cys', 'cysteine'),
        ('aro', 'aromatic (F, W, Y)'), ('hyd', 'aliphatic / hydrophobic'), ('pro', 'proline'), ('st', 'serine, threonine'),
        ('qn', 'glutamine, asparagine'), ('gly', 'glycine'), ('fcr', 'fraction charged')]
COLS = [('majority', 'majority'), ('half', 'half'), ('minority', 'minority')]   # all isoform-specific positions are in Fig. 5d
fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.9), gridspec_kw=dict(wspace=0.08))
lim = np.nanpercentile(np.abs(C.diff_pp), 97)
norm = TwoSlopeNorm(0, -lim, lim)
for ax, pos, t in zip(axes, ['internal', 'terminal'], 'ab'):
    M = np.full((len(FEAT), len(COLS)), np.nan); Q = np.ones_like(M); N = []
    for j, (k, _) in enumerate(COLS):
        z = C[C.kind == f'{k}|{pos}'].set_index('feature'); N.append(int(z.n_chunks.iloc[0]) if len(z) else 0)
        for i, (f, _) in enumerate(FEAT):
            if f in z.index: M[i, j] = z.loc[f, 'diff_pp']; Q[i, j] = z.loc[f, 'q_BH']
    im = ax.imshow(M, cmap='RdBu_r', norm=norm, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f'{M[i, j]:+.1f}', ha='center', va='center', fontsize=5.6, color='white' if abs(M[i, j]) > 0.6 * lim else INK)
                if Q[i, j] < 0.05: ax.plot(j + 0.38, i - 0.3, '.', ms=3, color=INK)
    ax.set_xticks(range(len(COLS))); ax.set_xticklabels([f'{l}\n({n:,})' for (_, l), n in zip(COLS, N)], fontsize=6)
    ax.axvline(2.5, color='white', lw=2)
    ax.set_yticks(range(len(FEAT))); ax.set_yticklabels([l for _, l in FEAT] if t == 'a' else [], fontsize=6.2)
    ax.set_title(f'{pos} isoform-specific vs shared disorder', fontsize=7)
    ax.set_xlabel("included in … of the gene's isoforms\n(disordered segments)"); tag(ax, t, -0.08 if t == 'b' else -0.55, 1.08)
    for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=axes, shrink=0.7, pad=0.02); cb.set_label('isoform-specific − shared (percentage points)', fontsize=6.5); cb.ax.tick_params(labelsize=6)
save(fig, 'extended_data/ed6_composition')
