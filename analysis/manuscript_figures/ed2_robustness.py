"""Extended Data Figure 2: robustness of the main results (disorder call, alignability, predictor agreement, estimators).

a, among internal exons whose splice sites are kept in at least 396 of 405 genomes, conservation of the 30 exonic and the 30
intronic bases next to the splice sites, by exon IDR fraction, relative to fully ordered exons: exonic conservation falls with
disorder, intronic conservation does not, so exonic conservation is a consequence of disorder rather than a measure of
alignability; b-d, three results recomputed under other disorder calls: alternative minus constitutive exon IDR fraction
(resampling genes), the within-gene effect of an exon with a poorly conserved splice site (gene fixed effects), and the change in canonical %IDR
per 100 Myr of gene age within families (family fixed effects); 95% intervals. Calls: the stored MetaPredict v3 IDRs;
MetaPredict v3 at disorder thresholds 0.4 and 0.6, minimum IDR length 10, 30 or 50, or raw scores without smoothing;
MetaPredict V1; MobiDB-lite consensus regions from UniProt (canonical isoforms only). Numbers from analysis/sensitivity/sensitivity.py and
analysis/exon_age/exon_age_robustness.py.
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from style import *

S = json.load(open(AN / 'sensitivity/outputs/sensitivity.json')); ER = json.load(open(AN / 'exon_age/outputs/robustness.json'))
CALLS = [('stored', 'stored (v3, paper)', BLUE), ('v3_dom_t0.4', 'v3, threshold 0.4', BLUE), ('v3_dom_t0.6', 'v3, threshold 0.6', BLUE),
         ('v3_dom_min10', 'v3, IDR ≥10 aa', BLUE), ('v3_dom_min30', 'v3, IDR ≥30 aa', BLUE), ('v3_dom_min50', 'v3, IDR ≥50 aa', BLUE),
         ('v3_raw_t0.5', 'v3, raw scores ≥0.5', BLUE), ('v1_dom_t0.3', 'MetaPredict V1', ORANGE), ('v1_raw_t0.3', 'V1, raw scores ≥0.3', ORANGE),
         ('mobidb_lite|canonical', 'MobiDB-lite (canonical)', ORANGE)]
CALLS = [c for c in CALLS if c[0] in S]
fig = plt.figure(figsize=(7.1, 7.6)); G = fig.add_gridspec(2, 1, height_ratios=[1.45, 1], hspace=0.42)
gs = G[0].subgridspec(1, 5, width_ratios=[0.85, 0.75, 1, 1, 1], wspace=0.4)
# a
ax = fig.add_subplot(gs[0, 0]); Y = ER['Y1_mean_exonic_flank_deep_by_idr']; ks = list(Y)
x = np.arange(len(ks))
for col, lab, c in [('eflank_weak', 'exonic bases', ORANGE), ('iflank_weak', 'intronic bases', BLUE)]:
    v = np.array([Y[k][col] for k in ks]); ax.plot(x, v - v[0], 'o-', color=c, ms=3.5, lw=1.2, label=lab)
refline(ax, 0, 'y', 'as in ordered exons', ha='left', x=0.02)
ax.set_xticks(x); ax.set_xticklabels(['0', '0–0.5', '0.5–1', '1'], fontsize=5.8); ax.set_xlabel('exon IDR fraction')
ax.set_ylabel('genomes keeping the flanking bases,\nrelative to fully ordered exons'); ax.legend(frameon=False, fontsize=5.6, loc='center left')
ax.set_title(f'{ER["Y1_deep_exons"]:,} exons with splice sites\nkept in ≥396 of 405 genomes', fontsize=6, color=MUTED, pad=4); grid(ax); tag(ax, 'a', -0.55)
# b-d forests
yy = np.arange(len(CALLS))[::-1]
for k, (key, i_est, xlab, ref) in enumerate([('S2_alt_minus_const', 0, 'alternative −\nconstitutive exons', 0),
                                              ('S3_young_exon_within_gene', 0, 'poorly conserved\nsplice site, within gene', 0),
                                              ('S4_slope_family_FE', 0, 'Δ %IDR per 100 Myr,\nwithin families', 0)]):
    ax = fig.add_subplot(gs[0, k + 2])
    for yi, (c, lab, col) in zip(yy, CALLS):
        v = S[c][key]; ax.plot([v[1], v[2]], [yi, yi], color=col, lw=1); ax.plot(v[0], yi, 'o', ms=3.6, color=col, mec='white', mew=0.4)
    refline(ax, ref, 'x', 'no effect' if k < 2 else 'no change', ha='left')
    ax.set_yticks(yy); ax.set_yticklabels([lab for _, lab, _ in CALLS] if k == 0 else [], fontsize=5.6); ax.set_xlabel(xlab); grid(ax, 'x')
    for t in ax.get_yticklabels():
        t.set_color(ORANGE if any(t.get_text() == l and col == ORANGE for _, l, col in CALLS) else INK)
    tag(ax, 'bcd'[k], -0.85 if k == 0 else -0.12)
ax.text(1.0, -0.2, 'blue: MetaPredict v3 settings · orange: other predictors', transform=ax.transAxes, ha='right', fontsize=5.4, color=ORANGE)
# e MetaPredict vs MobiDB-lite per class (formerly Fig. 3c); f estimators of the exon difference (formerly Fig. 4e)
r2 = G[1].subgridspec(1, 2, width_ratios=[0.9, 1.1], wspace=0.25)
PCK = json.load(open(AN / 'protein_classes/outputs/key_numbers.json')); CL = PCK['classes_by_median_pct_idr']; mb = PCK['K_mobidb_vs_metapredict']
ax = fig.add_subplot(r2[0])
for i, c in enumerate(CL):
    mp, md = mb['class_medians'][c]; ax.plot([i, i], [md, mp], color=LIGHT, lw=1.2); ax.plot(i, mp, 'o', ms=4, color=BLUE); ax.plot(i, md, 's', ms=3.4, mfc='white', mec=BLUE, mew=0.8)
ax.set_xticks(range(len(CL))); ax.set_xticklabels([SHORT[c] for c in CL], rotation=45, ha='right', fontsize=5.4); ax.set_ylabel('median % disordered'); grid(ax)
ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', color=BLUE, ms=4, label='MetaPredict IDRs'), plt.Line2D([], [], marker='s', ls='', mfc='white', mec=BLUE, ms=3.4, label='MobiDB-lite (UniProt)')],
          frameon=False, fontsize=5.8, loc='upper left')
ax.text(0.98, 0.04, f'per-protein ρ = {mb["spearman"]:.2f}', transform=ax.transAxes, ha='right', fontsize=5.8, color=MUTED); tag(ax, 'e', -0.2)
F3 = json.load(open(AN / 'manuscript_figures/outputs/fig3_numbers.json')); RF = json.load(open(AN / 'reference_free/outputs/key_numbers.json'))
est = [('raw difference, resampling genes', *F3['boot']), ('adjusted for length, position, transcripts', *F3['ols_adj'][:3]),
       ('within gene (fixed effects)', *F3['fixed_effects'][:3]), ('internal exons only, adjusted', *F3['internal_only']['adj']),
       ('usage on coding coordinates', F3['cds_level']['diff'], *F3['cds_level']['diff_ci95']),
       ('usage over all GENCODE transcripts, adjusted', *F3['gencode_wide_adj'][:3]),
       ('residues: isoform-specific − shared positions', *RF['pooled_diff_variable_minus_core'])]
ax = fig.add_subplot(r2[1]); yy = np.arange(len(est))[::-1]
for yi, (l, p, lo, hi) in zip(yy, est):
    ax.plot([lo, hi], [yi, yi], color=INK, lw=1); ax.plot(p, yi, 'o', ms=4, color=BLUE, mec='white', mew=0.5)
ax.set_yticks(yy); ax.set_yticklabels([e[0] for e in est], fontsize=5.8); ax.yaxis.tick_right(); ax.spines['right'].set_visible(True); ax.spines['left'].set_visible(False)
ax.set_xlim(0, 0.2); ax.set_xticks([0, 0.05, 0.1, 0.15, 0.2]); ax.set_xlabel('alternative − constitutive exon IDR fraction,\nall genes, stored call (95% CI)'); grid(ax, 'x'); tag(ax, 'f', -0.05)
save(fig, 'extended_data/ed2_robustness')
