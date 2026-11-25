"""Main Figure 5: how splicing changes IDRs between isoforms, with residues shared only when the same genomic codon encodes them.

a, fraction of coding positions in IDRs by inclusion (present in all, a majority, half or a minority of a gene's isoforms);
b, outcome of each IDR of one isoform in each other isoform of the gene (all ordered pairs), by IDR length;
c, outcomes by protein class;
d, composition of isoform-specific disorder minus the gene's shared disorder, in two panels: internal (flanked by shared
sequence) and terminal (beyond the first or last shared position; all, and in a minority of isoforms); rows grouped into
charge, linker/polar, hydrophobic and other; filled, Benjamini-Hochberg q < 0.05. The full grid is Extended Data Fig. 6.
"""
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from style import *

RF = json.load(open(AN / 'reference_free/outputs/key_numbers.json')); FM = json.load(open(AN / 'reference_free/outputs/fate_method_check.json'))
FR = json.load(open(AN / 'reference_free/outputs/fate_robustness.json'))
PC = json.load(open(AN / 'protein_classes/outputs/key_numbers.json')); CL = PC['classes_by_median_pct_idr']
H = pd.read_csv(AN / 'protein_classes/outputs/H_fates_by_class.csv', index_col=0)
R6 = pd.read_csv(AN / 'reanalysis_2026-10-26/outputs/R6_by_length.csv')
FATES = [('intact', LIGHT), ('internal deletion', BLUE), ('truncated', AQUA), ('removed', ORANGE)]
fig = plt.figure(figsize=(7.1, 7.0)); gs = fig.add_gridspec(2, 3, hspace=0.6, wspace=0.6, width_ratios=[0.8, 1.2, 0.9])
gs2 = gs[1, :].subgridspec(1, 3, width_ratios=[0.3, 1.1, 1.35], wspace=0.7)

def stacked(ax, rows, labels, fs=6.2):
    y = np.arange(len(rows))[::-1]; left = np.zeros(len(rows))
    for f, c in FATES:
        v = np.array([r.get(f, 0) for r in rows]); ax.barh(y, v, left=left, color=c, height=0.7, edgecolor='white', lw=0.5)
        for yi, l, vi in zip(y, left, v):
            if vi >= 2.5: ax.text(l + vi / 2, yi, f'{vi:.0f}', ha='center', va='center', fontsize=5.4 if vi >= 8 else (4.6 if vi >= 5 else 3.8),
                                  color=INK if f in ('intact', 'truncated') else 'white')   # every segment labelled; smaller type in narrow ones
        left += v
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=fs); ax.set_xlim(0, 100)

# a
ax = fig.add_subplot(gs[0, 0]); P = RF['pooled_frac_disordered_by_inclusion']; N = RF['pooled_positions_by_inclusion']
order = [('core', 'all'), ('majority', 'majority'), ('half', 'half'), ('minority', 'minority')]
for i, (k, l) in enumerate(order):
    ax.bar(i, P[k], color=MUTED if k == 'core' else ORANGE, alpha=1 if k in ('core', 'minority') else 0.7, width=0.7)
    ax.text(i, P[k] + 0.01, f'{100 * P[k]:.0f}%', ha='center', fontsize=6)
ax.set_xticks(range(4)); ax.set_xticklabels([l for _, l in order], fontsize=6, rotation=30, ha='right'); ax.set_xlabel("present in … of the gene's isoforms")
ax.set_ylabel('coding positions in IDRs (fraction)'); ax.set_ylim(0, 0.6); grid(ax); tag(ax, 'a', -0.45)
ax.text(0.02, 0.98, f'{RF["genes"]:,} genes\n{RF["positions"] / 1e6:.1f} M positions', transform=ax.transAxes, va='top', fontsize=5.6, color=MUTED)
# b
ax = fig.add_subplot(gs[0, 1:]); L = RF['all_pairs_fate_pct_by_length']
stacked(ax, [RF['all_pairs_fate_pct_pooled']] + [L[k] for k in ['20-50', '51-100', '101-200', '>200']], ['all IDRs', '20–50 aa', '51–100 aa', '101–200 aa', '>200 aa'])
ax.set_xlabel('IDRs of one isoform, followed\ninto each other isoform (%)'); tag(ax, 'b', -0.25)
fig.legend([Patch(color=c) for _, c in FATES], ['IDR in the other isoform: ' + FATES[0][0]] + [f for f, _ in FATES[1:]], ncol=4, frameon=False, fontsize=6.4, loc='lower center', bbox_to_anchor=(0.55, 0.93), handlelength=1, columnspacing=1.2)
# c
ax = fig.add_subplot(gs2[0, 1])
stacked(ax, [H.loc[c].to_dict() for c in CL[::-1]], [SHORT[c] for c in CL[::-1]], fs=5.6)
ax.set_xlabel('IDR outcomes by protein class (%)'); tag(ax, 'c', -0.55)
# d composition: internal and terminal isoform-specific disorder in separate panels, rows grouped by chemistry
C = pd.read_csv(AN / 'reference_free/outputs/composition_vs_core.csv')
GROUPS = [('charge', [('acid', 'acidic'), ('basic', 'basic'), ('ncpr', 'net charge')]),
          ('linker / polar', [('st', 'Ser, Thr'), ('pro', 'Pro'), ('qn', 'Gln, Asn'), ('gly', 'Gly')]),
          ('hydrophobic', [('hyd', 'aliphatic'), ('aro', 'aromatic')]),
          ('other', [('his', 'His'), ('cys', 'Cys')])]
ROWS, ypos, yl, gtop = [], [], 0, []
for g, fs in GROUPS:
    gtop.append((g, yl))
    for f in fs: ROWS.append(f); ypos.append(yl); yl -= 1
    yl -= 1.1
ypos = np.array(ypos); GREYD = '#9a9893'
PANELS = [('internal', [('variable|internal', 'o', INK, 0)]),
          ('terminal', [('variable|terminal', 'o', INK, 0.17), ('minority|terminal', 'D', GREYD, -0.17)])]
sub = gs2[0, 2].subgridspec(1, 2, wspace=0.08)
for pi, (title, sets) in enumerate(PANELS):
    ax = fig.add_subplot(sub[0, pi])
    for k, mk, col, off in sets:
        z = C[C.kind == k].set_index('feature')
        for yi, (f, _) in zip(ypos, ROWS):
            v = z.loc[f]; sig = v.q_BH < 0.05
            ax.plot([v.lo, v.hi], [yi + off] * 2, color=col, lw=0.7)
            ax.plot(v.diff_pp, yi + off, mk, ms=3.4 if mk == 'o' else 3.0, color=col if sig else 'white', mec=col, mew=0.8, zorder=3)
    ax.axvline(0, color=REF, lw=0.6, ls='--', zorder=0)
    for g, top in gtop[1:]: ax.axhline(top + 1.05, color=GRID, lw=0.6)
    ax.set_xlim(-4.6, 4.3); ax.set_ylim(ypos.min() - 0.7, 1.0); grid(ax, 'x'); ax.tick_params(axis='x', labelsize=5.6)
    ax.set_yticks(ypos); ax.set_yticklabels([l for _, l in ROWS] if pi == 0 else [], fontsize=5.8)
    if pi == 0:
        for g, top in gtop: ax.text(-0.62, top + 0.62, g, transform=ax.get_yaxis_transform(), fontsize=5.2, color=MUTED, style='italic', ha='left', va='center')
    # icon: shared (black) and isoform-specific (orange) sequence, and where the isoform-specific part sits
    T = ax.transAxes; segs = [(0.18, 0.28, INK), (0.30, 0.22, ORANGE), (0.54, 0.28, INK)] if title == 'internal' else [(0.18, 0.40, INK), (0.60, 0.22, ORANGE)]
    for x0, w, c in segs: ax.add_patch(plt.Rectangle((x0, 1.055), w, 0.035, transform=T, color=c, clip_on=False, lw=0))
    ax.text(0.5, 1.11, title, transform=T, ha='center', va='bottom', fontsize=6.2)
    if pi == 0: tag(ax, 'd', -0.75, 1.2)
fig.text(0.81, 0.05, 'isoform-specific − shared disorder (percentage points, 95% CI)', ha='center', fontsize=6.4)
leg = [plt.Line2D([], [], marker='o', ls='', color=INK, ms=3.4, label='all isoform-specific'),
       plt.Line2D([], [], marker='D', ls='', color=GREYD, ms=3, label='in a minority of isoforms'),
       plt.Line2D([], [], marker='o', ls='', mfc='white', mec=INK, ms=3.4, label='open: not significant (q ≥ 0.05)')]
fig.legend(handles=leg, frameon=False, fontsize=5.4, loc='lower center', bbox_to_anchor=(0.81, 0.005), ncol=3, handlelength=0.8, columnspacing=0.8)
save(fig, 'main/fig5_isoforms')
