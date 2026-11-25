"""Figure 3: alternative exons encode more predicted disorder than constitutive exons, across all genes.

Computes all-gene statistics (no TF split) from analysis/reanalysis_2026-10-26/outputs/unique_exons.csv and writes
outputs/fig3_numbers.json. Panels: a, cumulative distributions of exon IDR fraction; b, gene-paired means; c, effect
by gene disorder class; d, effect estimates under alternative estimators and usage definitions.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]; RE = REPO / 'analysis/reanalysis_2026-10-26/outputs'
FIG = REPO / 'manuscript/figures'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'pdf.fonttype': 42, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': 0.6, 'axes.edgecolor': '#52514e', 'xtick.color': '#52514e', 'ytick.color': '#52514e'})
BLUE, ORANGE, AQUA, MUTED, GRID, INK = '#2a78d6', '#eb6834', '#1baf7a', '#52514e', '#e6e5e0', '#0b0b0b'
g = pd.read_csv(RE / 'unique_exons.csv')
K = {}; rng = np.random.default_rng(20261026)
K['n_exons'] = len(g); K['n_alt'] = int(g.alt.sum()); K['n_const'] = int((1 - g.alt).sum()); K['n_genes'] = int(g.gene.nunique())
K['mean_alt'] = round(g[g.alt == 1].idr_frac.mean(), 4); K['mean_const'] = round(g[g.alt == 0].idr_frac.mean(), 4)
K['frac_idr_exon_alt'] = round(g[g.alt == 1].idr_exon.mean(), 4); K['frac_idr_exon_const'] = round(g[g.alt == 0].idr_exon.mean(), 4)
K['bimodal_frac_le0.1_or_ge0.9'] = round(float(((g.idr_frac <= 0.1) | (g.idr_frac >= 0.9)).mean()), 3)
idx = g.groupby('gene').indices; genes = np.array(list(idx))
def diff(rows):
    z = g.iloc[rows]; return z[z.alt == 1].idr_frac.mean() - z[z.alt == 0].idr_frac.mean()
bs = [diff(np.concatenate([idx[x] for x in rng.choice(genes, len(genes))])) for _ in range(1000)]
K['boot'] = [round(K['mean_alt'] - K['mean_const'], 4), round(np.percentile(bs, 2.5), 4), round(np.percentile(bs, 97.5), 4)]
pg = g.groupby(['gene', 'alt']).idr_frac.mean().unstack().dropna()
K['paired'] = dict(n=len(pg), frac_alt_higher=round(float((pg[1] > pg[0]).mean()), 3), frac_equal=round(float((pg[1] == pg[0]).mean()), 3),
                   median_delta=round(float((pg[1] - pg[0]).median()), 4), wilcoxon_p=float(f'{stats.wilcoxon(pg[1], pg[0]).pvalue:.3g}'))
cl = lambda m: m.fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(g.gene)[0]}, disp=0)
m = cl(smf.ols('idr_frac ~ alt + log_len + C(pos) + log_ntx', g))
K['ols_adj'] = [round(m.params.alt, 4), round(m.conf_int().loc['alt', 0], 4), round(m.conf_int().loc['alt', 1], 4), float(f'{m.pvalues.alt:.3g}')]
ml = cl(smf.logit('idr_exon ~ alt + log_len + C(pos) + log_ntx', g))
K['logit_OR'] = [round(np.exp(ml.params.alt), 3), round(np.exp(ml.conf_int().loc['alt', 0]), 3), round(np.exp(ml.conf_int().loc['alt', 1]), 3)]
w = g[g.gene.isin(pg.index)].copy()
for c in ['idr_frac', 'alt', 'log_len']: w[c + '_w'] = w[c] - w.groupby('gene')[c].transform('mean')
mf = smf.ols('idr_frac_w ~ alt_w + log_len_w - 1', w).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(w.gene)[0]})
K['fixed_effects'] = [round(mf.params.alt_w, 4), round(mf.conf_int().loc['alt_w', 0], 4), round(mf.conf_int().loc['alt_w', 1], 4), float(f'{mf.pvalues.alt_w:.3g}')]
# internal exons only (excludes alternative first/last exons from alternative promoters and polyadenylation)
gi = g[g.pos == 'internal']; mi = smf.ols('idr_frac ~ alt + log_len + log_ntx', gi).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(gi.gene)[0]})
K['internal_only'] = dict(n=len(gi), mean_alt=round(gi[gi.alt == 1].idr_frac.mean(), 4), mean_const=round(gi[gi.alt == 0].idr_frac.mean(), 4),
                          adj=[round(mi.params.alt, 4), round(mi.conf_int().loc['alt', 0], 4), round(mi.conf_int().loc['alt', 1], 4)])
CDS = json.load(open(RE / 'cds_level_usage.json')); RK = json.load(open(RE / 'key_numbers.json'))
BC = json.load(open(REPO / 'analysis/idr_landscape/outputs/splicing_by_class.json'))
K['cds_level'] = CDS; K['gencode_wide_adj'] = RK['R5_gencode_wide']['alt']; K['inclusion_slope'] = RK['R5_inclusion_slope']['incl']
json.dump(K, open(HERE / 'outputs/fig3_numbers.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in K.items() if k != 'cds_level'}, indent=1))

fig = plt.figure(figsize=(7.1, 4.6)); gs = fig.add_gridspec(2, 3, hspace=0.6, wspace=0.7, height_ratios=[1, 0.85])
axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1:]), fig.add_subplot(gs[0, 2])]
def tag(ax, t): ax.text(-0.38, 1.06, t, transform=ax.transAxes, fontweight='bold', fontsize=9)
ax = axs[0]
for a, c, lab in [(0, MUTED, 'constitutive'), (1, ORANGE, 'alternative')]:
    x = np.sort(g[g.alt == a].idr_frac.values); ax.plot(x, np.arange(1, len(x) + 1) / len(x), color=c, lw=1.5, label=f'{lab}\n(n = {len(x):,})')
ax.set_xlabel('exon IDR fraction'); ax.set_ylabel('cumulative fraction of exons'); ax.legend(frameon=False, fontsize=5.8, loc='lower right')
ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True); tag(ax, 'a')
ax = axs[1]
ax.hexbin(pg[0], pg[1], gridsize=25, cmap='Blues', bins='log', mincnt=1, lw=0.1, extent=(0, 1, 0, 1)); ax.plot([0, 1], [0, 1], color=MUTED, lw=0.6, ls='--')
ax.set_xlabel('constitutive exons\n(gene mean IDR fraction)'); ax.set_ylabel('alternative exons\n(gene mean IDR fraction)'); ax.set_aspect('equal')
ax.text(0.03, 0.97, f'{K["paired"]["n"]:,} genes', transform=ax.transAxes, va='top', fontsize=5.8); tag(ax, 'b')
ax = axs[2]
est = [('raw difference, resampling genes', K['boot'][0], K['boot'][1], K['boot'][2]),
       ('adjusted (length, position,\ntranscripts), clustered', *K['ols_adj'][:3]),
       ('within gene (fixed effects)', *K['fixed_effects'][:3]),
       ('internal exons only, adjusted', *K['internal_only']['adj']),
       ('usage on coding coordinates', CDS['diff'], *CDS['diff_ci95']),
       ('usage over all GENCODE\ncoding transcripts, adjusted', *K['gencode_wide_adj'][:3])]
y = np.arange(len(est))[::-1]
for yi, (lab, p, lo, hi) in zip(y, est):
    ax.plot([lo, hi], [yi, yi], color=INK, lw=1); ax.plot(p, yi, 'o', ms=4, color=BLUE, mec='white', mew=0.5)
ax.set_yticks(y); ax.set_yticklabels([e[0].replace(chr(10), ' ') for e in est], fontsize=6.3); ax.axvline(0, color=MUTED, lw=0.5)
ax.set_xlim(0, 0.22); ax.set_xlabel('alternative − constitutive exon IDR fraction (95% CI)'); ax.grid(axis='x', color=GRID, lw=0.5); ax.set_axisbelow(True)
ax.text(-0.62, 1.06, 'd', transform=ax.transAxes, fontweight='bold', fontsize=9)
ax = axs[3]; order = ['very low', 'intermediate', 'high']; cols = [AQUA, MUTED, ORANGE]
e = BC['exon_alt_minus_const']
for i, (c, col) in enumerate(zip(order, cols)):
    ax.bar(i, e[c]['diff'], color=col, width=0.62); ax.plot([i, i], [e[c]['lo'], e[c]['hi']], color=INK, lw=1)
ax.set_xticks(range(3)); ax.set_xticklabels([f'very low\n{e["very low"]["genes"]:,}', f'interm.\n{e["intermediate"]["genes"]:,}', f'high\n{e["high"]["genes"]:,}'], fontsize=6)
ax.set_ylabel('alternative − constitutive\nexon IDR fraction'); ax.set_xlabel('gene disorder class (genes)'); ax.grid(axis='y', color=GRID, lw=0.5); ax.set_axisbelow(True); tag(ax, 'c')
(FIG / 'supplementary').mkdir(parents=True, exist_ok=True)
for ext in ['pdf', 'png']: fig.savefig(FIG / 'supplementary' / f's2_exons_by_disorder_class.{ext}', dpi=300, bbox_inches='tight')
print('written', FIG / 'fig3_exons.pdf')
