"""Extended Data Figure 4: validation of the codon-level IDR fates.

a, the same canonical IDRs in alternative isoforms classified by difflib (autojunk off, as in earlier drafts; autojunk on;
blocks of at least 5 or 10 residues), by global Needleman-Wunsch alignment (all identical aligned residues; identical runs of
at least 5) and by codon keys (first base, as in the paper; all three bases); b, fates over all ordered pairs pooled, with
genes weighted equally, without the 5% of genes contributing most comparisons, with one random canonical-alternative pair
per gene, and by the number of isoforms in the gene; c, by direction; d, fraction of an IDR's kept codons that are
disordered in the other isoform's own prediction, by fate. Numbers from analysis/reference_free/fate_robustness.py.
"""
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from style import *

FR = json.load(open(AN / 'reference_free/outputs/fate_robustness.json'))
R = pd.read_csv(AN / 'reference_free/outputs/fate_robustness_all_pairs.csv.gz', usecols=['fate', 'frac_disordered_in_b'])
FATES = [('intact', LIGHT), ('internal deletion', BLUE), ('truncated', AQUA), ('removed', ORANGE)]
fig = plt.figure(figsize=(7.1, 6.4)); gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.95, width_ratios=[1.1, 1])


def stacked(ax, rows, labels, fs=5.8):
    y = np.arange(len(rows))[::-1]; left = np.zeros(len(rows))
    for f, c in FATES:
        v = np.array([r.get(f, 0) for r in rows]); ax.barh(y, v, left=left, color=c, height=0.7, edgecolor='white', lw=0.5)
        for yi, l, vi in zip(y, left, v):
            if vi >= 2.5: ax.text(l + vi / 2, yi, f'{vi:.0f}', ha='center', va='center', fontsize=5.4 if vi >= 8 else (4.6 if vi >= 5 else 3.8),
                                  color=INK if f in ('intact', 'truncated') else 'white')   # every segment labelled; smaller type in narrow ones
        left += v
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=fs); ax.set_xlim(0, 100)


ax = fig.add_subplot(gs[0, 0]); P = FR['R1_fate_pct']; A = FR['R1_agreement_with_codon_full']
M = [('difflib', 'difflib, autojunk off'), ('difflib_autojunk', 'difflib, autojunk on'), ('difflib_min5', 'difflib, blocks ≥5'), ('difflib_min10', 'difflib, blocks ≥10'),
     ('needle', 'global alignment'), ('needle_min5', 'global, runs ≥5'), ('codon_first_base', 'genomic codon, first base'), ('codon_full', 'genomic codon, three bases')]
stacked(ax, [P[k] for k, _ in M], [l for _, l in M])
for yi, (k, _) in zip(np.arange(len(M))[::-1], M):
    ax.text(101, yi, f'{100 * A[k]:.0f}%', va='center', fontsize=5.2, color=MUTED)
ax.text(101, len(M) - 0.45, 'agree with\nthree-base match', fontsize=5, color=MUTED, va='bottom')
ax.set_xlabel(f'the same {FR["R1_idrs"]["with_needle"]:,} canonical IDRs in\nalternative isoforms (%)'); tag(ax, 'a', -0.62)
ax = fig.add_subplot(gs[0, 1]); S = FR['R3_by_isoform_count']
rows = [('pooled (paper)', FR['R3_pooled']), ('genes weighted equally', FR['R3_gene_weighted']), ('without top 5% of genes', FR['R3_without_top_5pct_genes']),
        ('one pair per gene', FR['R3_one_random_canonical_to_alternative_pair_per_gene'])] + \
       [(f'{k} isoforms ({S[k]["genes"]:,} genes)', S[k]) for k in ['2', '3', '4-5', '6+']]
stacked(ax, [r for _, r in rows], [l for l, _ in rows])
ax.set_xlabel('IDR comparisons (%)'); tag(ax, 'b', -0.7)
ax.text(0.0, 1.02, f'top 5% of genes give {100 * FR["R3_share_of_comparisons_top_5pct_genes"]:.0f}% of comparisons', transform=ax.transAxes, fontsize=5.4, color=MUTED)
ax = fig.add_subplot(gs[1, 0]); D = FR['R3_by_direction']
stacked(ax, [D['canonical_to_alternative'], D['alternative_to_canonical'], D['alternative_to_alternative']],
        ['canonical → alternative', 'alternative → canonical', 'alternative → alternative'])
ax.set_xlabel('IDR followed from one isoform into another (%)'); tag(ax, 'c', -0.62)
ax = fig.add_subplot(gs[1, 1]); bins = np.linspace(0, 1, 21)
for f, c in [('intact', MUTED), ('internal deletion', BLUE), ('truncated', AQUA)]:
    v = R[R.fate == f].frac_disordered_in_b.dropna()
    ax.hist(v, bins=bins, weights=np.full(len(v), 100 / len(v)), histtype='step', color=c, lw=1.2, label=f'{f} ({100 * (v >= 0.9).mean():.0f}% ≥0.9)')
ax.set_yscale('log'); ax.set_xlabel("kept codons disordered in the\nother isoform's prediction (fraction)"); ax.set_ylabel('IDR comparisons (%)')
ax.legend(frameon=False, fontsize=5.4, loc='upper left'); grid(ax); tag(ax, 'd', -0.3)
fig.legend([Patch(color=c) for _, c in FATES], ['IDR in the other isoform: ' + FATES[0][0]] + [f for f, _ in FATES[1:]], ncol=4, frameon=False, fontsize=6.2,
           loc='lower center', bbox_to_anchor=(0.5, 0.95), handlelength=1, columnspacing=1.2)
save(fig, 'extended_data/ed4_fate_validation')
