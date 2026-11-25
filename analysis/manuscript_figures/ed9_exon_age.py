"""Extended Data Figure 9: exons with poorly conserved splice sites encode more disorder, within genes.

Exon conservation = the lower of its two splice sites' counts of the 405 aligned mammalian genomes that keep the GT/AG
dinucleotide (Minkin & Salzberg 2025); internal coding exons only. a, genes with both poorly conserved (<100) and deeply conserved (>=350) internal exons; b, the within-gene
effect of an exon with a poorly conserved splice site, re-estimated to answer simple alternative explanations; d, composition of disordered
residues from poorly conserved minus deeply conserved exons of the same genes, next to young-minus-old genes (core disorder).
"""
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from style import *

O = AN / 'exon_age/outputs'; K = json.load(open(O / 'key_numbers.json'))
D = pd.read_csv(O / 'exon_conservation.csv'); I = D[(D.pos == 'internal') & (D.n_sites == 2)]
EV = json.load(open(AN / 'evolution/outputs/key_numbers.json'))['E3_young_minus_old_core_pp']
fig = plt.figure(figsize=(7.1, 3.3)); gs = fig.add_gridspec(1, 2, wspace=1.05, width_ratios=[0.8, 1.2], left=0.1, right=0.98, top=0.88, bottom=0.18)
I = I.copy()
# b
ax = fig.add_subplot(gs[0, 0]); P = K['X2_paired_genes']
pg = I[I.cons.notna()].assign(b=np.where(I.cons < 100, 'young', np.where(I.cons >= 350, 'old', None))).dropna(subset=['b'])
pg = pg.groupby(['gene', 'b']).idr_frac.mean().unstack().dropna()
ax.hexbin(pg.old, pg.young, gridsize=18, cmap='Blues', bins='log', mincnt=1, extent=(0, 1, 0, 1), lw=0.1)
ax.plot([0, 1], [0, 1], '--', color=REF, lw=0.7); ax.text(0.62, 0.66, 'equal', rotation=45, fontsize=5.4, color=MUTED, ha='center'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xlabel('deeply conserved exons (≥350)\nmean IDR fraction'); ax.set_ylabel('poorly conserved exons (<100)\nmean IDR fraction')
ax.text(0.97, 0.05, f'{P["genes"]:,} genes\n{100 * P["frac_young_higher"]:.0f}% above, {100 * P["frac_young_lower"]:.0f}% below\nWilcoxon P = {P["wilcoxon_p"]:.0e}',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=5.6); tag(ax, 'a', -0.3)
# c: the within-gene result, asked as a series of plain questions
ax = fig.add_subplot(gs[0, 1]); ER = json.load(open(AN / 'exon_age/outputs/robustness.json')); RM = {m['model']: m for m in ER['models']}
S = json.load(open(AN / 'sensitivity/outputs/sensitivity.json')); Y4 = ER['Y4_matched']
f = lambda m: (RM[m]['coef'], RM[m]['lo'], RM[m]['hi'])
QS = [('Main result', [('all internal exons', (K['X2_fe_young']['coef'], K['X2_fe_young']['lo'], K['X2_fe_young']['hi']))]),
      ('Only because alternative exons are disordered?', [('constitutive exons only', K['X2_fe_young_constitutive_only'][:3]), ('alternative exons only', K['X2_fe_young_alternative_only'][:3])]),
      ('Only because fast-evolving sequence aligns poorly?', [('accounting for how well the nearby intron aligns', f('+ intronic flank covariate')), ('only exons whose nearby intron aligns well', f('intronic flank >= 200'))]),
      ('Only because of exon length or usage?', [('pairs of exons matched for length and usage', (Y4['diff'], Y4['lo'], Y4['hi']))]),
      ('Only because of the disorder predictor?', [('MetaPredict V1', S['v1_dom_t0.3']['S3_young_exon_within_gene'][:3]), ('MobiDB-lite (canonical exons)', S['mobidb_lite|canonical']['S3_young_exon_within_gene'][:3])])]
y, ticks, labs = 0, [], []
for q, rows in QS:
    ax.text(-0.02, y, q, transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=5.6, fontweight='bold')
    y -= 1.0
    for l, (p_, lo, hi) in rows:
        ax.plot([lo, hi], [y, y], color=INK, lw=1); ax.plot(p_, y, 'o', ms=4, color=BLUE, mec='white', mew=0.5); ticks.append(y); labs.append(l); y -= 0.85
    y -= 0.6
refline(ax, 0, 'x', 'no difference', ha='left'); ax.set_yticks(ticks); ax.set_yticklabels(labs, fontsize=5.6); ax.set_ylim(y + 0.9, 0.6); ax.set_xlim(-0.02, 0.16); ax.set_xticks([0, 0.05, 0.1, 0.15])
ax.set_xlabel('extra IDR fraction of a poorly conserved exon,\ncompared with the same gene’s other exons (95% CI)'); grid(ax, 'x'); tag(ax, 'b', -1.0, 1.08)
save(fig, 'extended_data/ed9_exon_age')
