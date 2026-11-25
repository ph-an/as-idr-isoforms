"""Robustness of the exon-age result (Fig. 6f-h; Extended Data Figs 4-5).

The headline estimate (exon_age.py, X2) is that, within genes, an internal exon whose less conserved splice site keeps its
GT/AG in fewer than 100 of 405 mammals encodes +0.107 more IDR fraction. That estimate shrinks when exons are restricted by
the conservation of their flanking EXONIC bases and reverses when exonic conservation is a covariate. This script asks
which of these controls is appropriate and whether the result survives other definitions.

  Y1  Is exonic conservation a confounder or a consequence of disorder? Among exons whose splice sites are conserved in
      >= 396 genomes (present across mammals by any definition), exonic-flank conservation should not depend on disorder
      if it only measured alignability. We compare it with the INTRONIC flank (the 30 intronic bases next to each site),
      which reflects whether the locus aligns but is not constrained by the protein.
  Y2  Within-gene models with the intronic flank as the alignability covariate, and restricted by intronic-flank
      thresholds, next to the original exonic-flank versions.
  Y3  Alternative definitions of a lineage-restricted exon: thresholds of 50, 100, 150 and 200 genomes; the acceptor or the
      donor alone; the more conserved site instead of the less conserved; a continuous count; and a "site deficit"
      definition, the splice site conserved in at least 150 fewer genomes than its own intronic flank (the locus aligns
      but the splice signal is absent).
  Y4  Young vs old exons matched within gene on length and usage (nearest old exon of the same gene and usage), paired
      differences with gene bootstrap.

Writes outputs/robustness.json and outputs/robustness_models.csv.
"""
import os, json
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); OUT = os.path.join(HERE, 'outputs')
rng = np.random.default_rng(20261105); K = {}
U = pd.read_parquet(os.path.join(OUT, 'sites_unique.parquet'))
# per-position columns run in transcript orientation: acceptor exon = positions 2..31 (flank_pos), intron = -30..-1 (flank_neg);
# donor exon = -30..-1 (flank_neg), intron = 2..31 (flank_pos)
U['intron_flank'] = np.where(U.site_type == 'a', U.flank_neg, U.flank_pos)
U['exon_flank'] = np.where(U.site_type == 'a', U.flank_pos, U.flank_neg)
ik = U.set_index(['chr', 'strand', 'pos', 'site_type']).intron_flank
ek = U.set_index(['chr', 'strand', 'pos', 'site_type']).exon_flank
E = pd.read_parquet(os.path.join(ROOT, 'data/genomic/processed/exon_level_table.parquet'), columns=['base_accession', 'exon_id', 'chr', 'start', 'end', 'strand', 'is_coding_exon'])
X = E[E.is_coding_exon].drop_duplicates(['base_accession', 'exon_id']).copy(); p = (X.strand == '+').values
acc = np.where(p, X.start - 2, X.end + 2); don = np.where(p, X.end + 1, X.start - 1)
lk = lambda s, pos, t: s.reindex(pd.MultiIndex.from_arrays([X.chr.values, X.strand.values, pos, np.array([t] * len(X))])).values
X['ifl_acc'], X['ifl_don'] = lk(ik, acc, 'a'), lk(ik, don, 'd'); X['efl_acc'], X['efl_don'] = lk(ek, acc, 'a'), lk(ek, don, 'd')
D = pd.read_csv(os.path.join(OUT, 'exon_conservation.csv'))
D = D.merge(X[['base_accession', 'exon_id', 'ifl_acc', 'ifl_don', 'efl_acc', 'efl_don']], on=['base_accession', 'exon_id'], how='left')
I = D[(D.pos == 'internal') & (D.n_sites == 2)].copy()
# flank of the less conserved site (the site that defines the exon's count), and the min over both sites
weak_is_acc = I.cons_acc <= I.cons_don
I['iflank'] = np.fmin(I.ifl_acc, I.ifl_don); I['iflank_weak'] = np.where(weak_is_acc, I.ifl_acc, I.ifl_don)
I['eflank_weak'] = np.where(weak_is_acc, I.efl_acc, I.efl_don)
I['iflank100'] = I.iflank / 100; I['flank100'] = I.flank / 100
K['exons'] = len(I); K['genes'] = int(I.gene.nunique())
K['flank_medians_by_bin'] = I.groupby('cbin')[['flank', 'iflank']].median().round(1).to_dict(orient='index')


def fe(df, y, xs):
    z = df.dropna(subset=[y] + xs).copy(); z = z[z.groupby('gene')[y].transform('size') > 1]
    for c in [y] + xs: z[c + '_d'] = z[c] - z.groupby('gene')[c].transform('mean')
    m = smf.ols(f'{y}_d ~ ' + ' + '.join(c + '_d' for c in xs) + ' - 1', z).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(z.gene)[0]})
    return m, z


rows = []
def fit(label, df, var='young', extra=()):
    xs = [var, 'log_len', 'alt', 'log_ntx'] + list(extra)
    m, z = fe(df, 'idr_frac', xs); ci = m.conf_int().loc[var + '_d']
    r = dict(model=label, var=var, coef=round(m.params[var + '_d'], 4), lo=round(ci[0], 4), hi=round(ci[1], 4), p=float(f'{m.pvalues[var + "_d"]:.3g}'),
             genes=int(z.gene.nunique()), exons=len(z), n_young=int(z[var].sum()) if set(z[var].unique()) <= {0, 1} else None)
    rows.append(r); return r

# ---- Y1
deep = I[I.cons >= 396]
K['Y1_deep_exons'] = len(deep)
K['Y1_spearman_idr_vs_exonic_flank_deep'] = [round(v, 4) for v in stats.spearmanr(deep.idr_frac, deep.eflank_weak, nan_policy='omit')]
K['Y1_spearman_idr_vs_intronic_flank_deep'] = [round(v, 4) for v in stats.spearmanr(deep.idr_frac, deep.iflank_weak, nan_policy='omit')]
K['Y1_mean_exonic_flank_deep_by_idr'] = deep.groupby(pd.cut(deep.idr_frac, [-0.01, 0.0, 0.5, 0.99, 1.0], labels=['0', '(0,0.5]', '(0.5,1)', '1']), observed=True)[['eflank_weak', 'iflank_weak']].mean().round(1).to_dict(orient='index')
m1, _ = fe(deep.assign(ef=deep.eflank_weak / 100, inf=deep.iflank_weak / 100), 'ef', ['idr_frac', 'log_len', 'alt'])
K['Y1_within_gene_exonic_flank_per_unit_idr_deep'] = [round(m1.params.idr_frac_d * 100, 2), float(f'{m1.pvalues.idr_frac_d:.3g}')]
m1i, _ = fe(deep.assign(inf=deep.iflank_weak / 100), 'inf', ['idr_frac', 'log_len', 'alt'])
K['Y1_within_gene_intronic_flank_per_unit_idr_deep'] = [round(m1i.params.idr_frac_d * 100, 2), float(f'{m1i.pvalues.idr_frac_d:.3g}')]

# ---- Y2
fit('original (young <100)', I)
fit('+ exonic flank covariate', I, extra=['flank100'])
fit('+ intronic flank covariate', I, extra=['iflank100'])
for t in [100, 150, 200, 250]:
    fit(f'intronic flank >= {t}', I[I.iflank >= t])
for t in [200, 300, 350]:
    fit(f'exonic flank >= {t}', I[I.flank >= t])
# ---- Y3
for t in [50, 150, 200]:
    fit(f'young < {t}', I.assign(y2=(I.cons < t).astype(float)), 'y2')
fit('acceptor only < 100', I.assign(y2=(I.cons_acc < 100).astype(float)), 'y2')
fit('donor only < 100', I.assign(y2=(I.cons_don < 100).astype(float)), 'y2')
fit('more conserved site < 100', I.assign(y2=(I[['cons_acc', 'cons_don']].max(axis=1) < 100).astype(float)), 'y2')
fit('continuous: per 100 fewer genomes', I.assign(y2=(405 - I.cons) / 100), 'y2')
fit('continuous + intronic flank', I.assign(y2=(405 - I.cons) / 100), 'y2', ['iflank100'])
I['deficit'] = I.iflank_weak - I.cons
fit('site deficit >= 150 (locus aligns, site absent)', I.assign(y2=(I.deficit >= 150).astype(float)), 'y2')
fit('site deficit >= 150, intronic flank >= 200', I[I.iflank_weak >= 200].assign(y2=lambda z: (z.deficit >= 150).astype(float)), 'y2')
K['Y3_site_deficit_exons'] = int((I.deficit >= 150).sum())

# ---- Y4 matched pairs within gene
pairs = []
for g, z in I.dropna(subset=['young']).groupby('gene'):
    yv, ov = z[z.young == 1], z[z.young == 0]
    if yv.empty or ov.empty: continue
    for r in yv.itertuples():
        cand = ov[ov.alt == r.alt]
        if cand.empty: continue
        j = (cand.log_len - r.log_len).abs().idxmin()
        if abs(cand.log_len[j] - r.log_len) <= np.log10(1.5): pairs.append((g, r.idr_frac, cand.idr_frac[j], r.alt))
P = pd.DataFrame(pairs, columns=['gene', 'young', 'old', 'alt'])
gm = P.groupby('gene')[['young', 'old']].mean(); dd = gm.young - gm.old
bs = [dd.sample(len(dd), replace=True, random_state=int(s)).mean() for s in rng.integers(0, 2 ** 31, 1000)]
K['Y4_matched'] = dict(pairs=len(P), genes=len(gm), mean_young=round(gm.young.mean(), 4), mean_old=round(gm.old.mean(), 4), diff=round(dd.mean(), 4),
                       lo=round(np.percentile(bs, 2.5), 4), hi=round(np.percentile(bs, 97.5), 4), wilcoxon_p=float(f'{stats.wilcoxon(gm.young, gm.old).pvalue:.3g}'))
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, 'robustness_models.csv'), index=False); K['models'] = rows
json.dump(K, open(os.path.join(OUT, 'robustness.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps({k: v for k, v in K.items() if k != 'models'}, indent=1, default=str)); print(R.to_string())
