"""Evolutionary analyses of intrinsic disorder and its alternative splicing (follow-up of Song & Akey 2026).

Gene age is the GenOrigin estimate (Ma; bins <100, 100-500, 500-1000, >1000). Protein classes come from
analysis/protein_classes; core/variable positions and IDR fates from analysis/reference_free (shared genomic codons).

  E1  disorder-to-order by protein class: canonical %IDR vs gene age (Spearman per class); OLS of gene age on %IDR with
      a class interaction, adjusted for log length (the form of Song & Akey's model, with classes instead of TF/non-TF)
  E2  splicing across gene age: isoforms per gene, alternative-exon fraction, alternative - constitutive exon IDR fraction,
      variable - core disorder, and codon-level IDR fates, by age bin (overall and per class)
  E3  IDR composition across gene age: composition of core disorder by age bin (positive, aromatic, acidic ...), and the
      composition of variable disorder (by inclusion) relative to core disorder within each age bin; similarity between
      the 'young vs old' and 'rarely included vs core' composition contrasts
  E4  disease inheritance: canonical %IDR, splicing coupling and variable - core disorder in genes with autosomal
      dominant vs autosomal recessive disease annotations (UniProt 'Involvement in disease')
Writes outputs/key_numbers.json and CSV tables.
"""
import os, json
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats

ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
HERE = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(HERE, '..'); OUT = os.path.join(HERE, 'outputs'); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261104)
AGE = ['<100', '100-500', '500-1000', '>1000']
G = pd.read_csv(os.path.join(A, 'protein_classes/outputs/gene_annotations.csv'), index_col=0)
CLASSES = G.groupby('pclass').pct_idr.median().sort_values().index.tolist()
K = {}

# ---- E1 disorder vs age by class
ga = G.dropna(subset=['gene_age']).copy(); ga['log_len'] = np.log10(ga.length)
K['E1_genes_with_age'] = len(ga)
rows = []
for c in ['All'] + CLASSES:
    d = ga if c == 'All' else ga[ga.pclass == c]
    r = stats.spearmanr(d.gene_age, d.pct_idr)
    med = d.groupby('age_bin').pct_idr.median().reindex(AGE)
    kw = stats.kruskal(*[g.pct_idr for _, g in d.groupby('age_bin') if len(g) > 5])
    rows.append(dict(pclass=c, genes=len(d), rho=r[0], p=r[1], kruskal_p=kw.pvalue, **{f'median_{a}': med[a] for a in AGE},
                     **{f'n_{a}': int((d.age_bin == a).sum()) for a in AGE}))
T1 = pd.DataFrame(rows); T1.to_csv(os.path.join(OUT, 'E1_disorder_vs_age_by_class.csv'), index=False)
K['E1'] = T1.round(4).to_dict(orient='records')
# TF labels (HumanTFs) as in Song & Akey
for lab, d in [('HumanTFs TF', ga[ga.is_tf == True]), ('non-TF', ga[ga.is_tf != True])]:
    r = stats.spearmanr(d.gene_age, d.pct_idr); K[f'E1_{lab}'] = [round(r[0], 3), float(f'{r[1]:.3g}'), len(d)]
m = smf.ols("gene_age ~ pct_idr * C(pclass, Treatment('Metabolic enzymes')) + log_len", ga).fit(cov_type='HC3')
sl = {}
for c in CLASSES:
    if c == 'Metabolic enzymes':
        b = m.params['pct_idr']; se = m.bse['pct_idr']
    else:
        t = f"pct_idr:C(pclass, Treatment('Metabolic enzymes'))[T.{c}]"
        b = m.params['pct_idr'] + m.params[t]; cov = m.cov_params(); se = np.sqrt(cov.loc['pct_idr', 'pct_idr'] + cov.loc[t, t] + 2 * cov.loc['pct_idr', t])
    sl[c] = [round(b, 3), round(b - 1.96 * se, 3), round(b + 1.96 * se, 3)]
K['E1_age_on_pct_idr_slope_by_class_Ma_per_pct'] = sl

# ---- E2 splicing across age
ux = pd.read_csv(os.path.join(A, 'reanalysis_2026-10-26/outputs/unique_exons.csv')); ux['age_bin'] = ux.gene.map(G.age_bin); ux['pclass'] = ux.gene.map(G.pclass)
uidx = ux.groupby('gene').indices
def ediff(z): return z[z.alt == 1].idr_frac.mean() - z[z.alt == 0].idr_frac.mean()
def boot(genes, B=500):
    pt = ediff(ux.iloc[np.concatenate([uidx[g] for g in genes])])
    bs = [ediff(ux.iloc[np.concatenate([uidx[g] for g in rng.choice(genes, len(genes))])]) for _ in range(B)]
    return pt, np.percentile(bs, 2.5), np.percentile(bs, 97.5)
RF = pd.read_csv(os.path.join(A, 'reference_free/outputs/gene_core_variable.csv')).set_index('gene')
RF['age_bin'] = RF.index.map(G.age_bin); RF['pclass'] = RF.index.map(G.pclass)
F = pd.read_csv(os.path.join(A, 'reference_free/outputs/all_pairs_idr_fates.csv')); F['age_bin'] = F.gene.map(G.age_bin)
afr = ux.groupby('gene').alt.mean()
rows = []
for a in AGE:
    gs = ux[ux.age_bin == a].gene.unique(); d, lo, hi = boot(gs)
    rf = RF[RF.age_bin == a].dropna(subset=['dis_core', 'dis_var'])
    fa = F[F.age_bin == a].groupby('fate').n.sum(); fa = 100 * fa / fa.sum()
    gg = G[G.age_bin == a]
    rows.append(dict(age_bin=a, genes=len(gg), mean_iso=gg.n_iso.mean(), frac_ge3_iso=(gg.n_iso >= 3).mean(), median_idr_range=gg.idr_range.median(),
                     alt_exon_fraction=afr[afr.index.isin(gs)].mean(), coupling=d, coupling_lo=lo, coupling_hi=hi,
                     dis_core=rf.dis_core.mean(), dis_var=rf.dis_var.mean(), var_minus_core=(rf.dis_var - rf.dis_core).mean(),
                     **{f'fate_{k}': fa.get(k, np.nan) for k in ['intact', 'internal deletion', 'truncated', 'removed']}))
T2 = pd.DataFrame(rows); T2.to_csv(os.path.join(OUT, 'E2_splicing_by_age.csv'), index=False); K['E2'] = T2.round(4).to_dict(orient='records')
K['E2_spearman_age_vs_n_iso'] = [round(x, 4) for x in stats.spearmanr(ga.gene_age, ga.n_iso)]
rows = []
for c in CLASSES:
    for a in AGE:
        gs = ux[(ux.pclass == c) & (ux.age_bin == a)].gene.unique()
        rf = RF[(RF.pclass == c) & (RF.age_bin == a)].dropna(subset=['dis_core', 'dis_var'])
        if len(gs) < 25: continue
        d, lo, hi = boot(gs, B=300)
        rows.append(dict(pclass=c, age_bin=a, genes=len(gs), coupling=d, lo=lo, hi=hi, var_minus_core=(rf.dis_var - rf.dis_core).mean(),
                         median_pct_idr=G[(G.pclass == c) & (G.age_bin == a)].pct_idr.median()))
T2c = pd.DataFrame(rows); T2c.to_csv(os.path.join(OUT, 'E2_splicing_by_class_age.csv'), index=False)

# ---- E3 composition across age
Ch = pd.read_csv(os.path.join(A, 'reference_free/outputs/disorder_chunks_by_inclusion.csv')); Ch['age_bin'] = Ch.gene.map(G.age_bin)
Ch['pclass'] = Ch.gene.map(G.pclass); Ch['pos'] = Ch.basic + Ch.his      # K, R, H as in Song & Akey
FE = ['aro', 'acid', 'basic', 'his', 'pos', 'pro', 'st', 'qn', 'gly', 'hyd', 'cys', 'ncpr']
core = Ch[Ch.kind == 'core'].dropna(subset=['age_bin'])
T3 = core.groupby('age_bin')[FE].mean().reindex(AGE); T3['n'] = core.groupby('age_bin').size().reindex(AGE)
T3.to_csv(os.path.join(OUT, 'E3_core_disorder_composition_by_age.csv')); K['E3_core_composition_by_age'] = (T3 * 100).round(2).to_dict(orient='index')
K['E3_kruskal_core_by_age'] = {f: float(f'{stats.kruskal(*[g[f] for _, g in core.groupby("age_bin")]).pvalue:.3g}') for f in FE}
# young-vs-old contrast of core disorder (gene-clustered OLS, log size adjusted)
core = core.assign(young=(core.age_bin.isin(['<100', '100-500'])).astype(int), old=(core.age_bin == '>1000').astype(int), logn=np.log10(core.n))
yo = core[(core.young == 1) | (core.old == 1)]
young_vs_old = {}
for f in FE:
    mm = smf.ols(f'{f} ~ young + logn', yo).fit(cov_type='HC3'); young_vs_old[f] = 100 * mm.params.young
K['E3_young_minus_old_core_pp'] = {k: round(v, 3) for k, v in young_vs_old.items()}
rf_c = pd.read_csv(os.path.join(A, 'reference_free/outputs/composition_vs_core.csv'))
minority = rf_c[rf_c.kind == 'minority|terminal'].set_index('feature').diff_pp
variable = rf_c[rf_c.kind.isin(['variable|internal', 'variable|terminal'])].groupby('feature').diff_pp.mean()
common = [f for f in FE if f in minority.index]
yv = np.array([young_vs_old[f] for f in common])
K['E3_similarity'] = dict(features=common,
                          spearman_young_old_vs_minority_terminal=[round(x, 3) for x in stats.spearmanr(yv, minority.reindex(common))],
                          spearman_young_old_vs_variable_mean=[round(x, 3) for x in stats.spearmanr(yv, variable.reindex(common))])
# variable vs core composition within each age bin
rows = []
for a in AGE:
    x = Ch[(Ch.age_bin == a) & Ch.kind.isin(['core', 'variable|internal', 'variable|terminal'])].copy()
    x['s'] = (x.kind != 'core').astype(int); x['logn'] = np.log10(x.n)
    for f in ['aro', 'acid', 'pos', 'his', 'cys', 'ncpr', 'hyd', 'pro']:
        mm = smf.ols(f'{f} ~ s + logn', x).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(x.gene)[0]})
        rows.append(dict(age_bin=a, feature=f, diff_pp=100 * mm.params.s, lo=100 * mm.conf_int().loc['s', 0], hi=100 * mm.conf_int().loc['s', 1], p=mm.pvalues.s,
                         n_var=int(x.s.sum())))
T3b = pd.DataFrame(rows); T3b.to_csv(os.path.join(OUT, 'E3_variable_vs_core_by_age.csv'), index=False)

# ---- E4 disease inheritance
u = pd.read_csv(os.path.join(ROOT, 'data/raw/uniprot/uniprotkb_annotated'), sep='\t', usecols=['Entry', 'Involvement in disease'], low_memory=False).set_index('Entry')
txt = u['Involvement in disease'].reindex(G.index).fillna('').str.lower()
G['AD'] = txt.str.contains('autosomal dominant'); G['AR'] = txt.str.contains('autosomal recessive')
has_dz = u['Involvement in disease'].reindex(G.index).fillna('').str.strip() != ''
grp = pd.Series(np.select([G.AD & ~G.AR, G.AR & ~G.AD, G.AD & G.AR, ~has_dz], ['AD only', 'AR only', 'AD and AR', 'no disease text'], 'disease, inheritance not stated'), index=G.index)
rows = []
for lab in ['AD only', 'AR only', 'AD and AR', 'disease, inheritance not stated', 'no disease text']:
    gs = grp.index[grp == lab]; gx = [g for g in gs if g in uidx]
    d, lo, hi = boot(np.array(gx), B=300) if len(gx) > 30 else (np.nan,) * 3
    rf = RF.reindex(gs).dropna(subset=['dis_core', 'dis_var'])
    pv = G.loc[gs, 'pct_idr'].dropna().values; vc = (rf.dis_var - rf.dis_core).values      # gene-bootstrap intervals for the other two measures
    bm = [np.median(rng.choice(pv, len(pv))) for _ in range(1000)]; bv = [rng.choice(vc, len(vc)).mean() for _ in range(1000)]
    rows.append(dict(group=lab, genes=len(gs), median_pct_idr=G.loc[gs, 'pct_idr'].median(), median_lo=np.percentile(bm, 2.5), median_hi=np.percentile(bm, 97.5),
                     mean_iso=G.loc[gs, 'n_iso'].mean(), coupling=d, lo=lo, hi=hi,
                     var_minus_core=vc.mean(), vc_lo=np.percentile(bv, 2.5), vc_hi=np.percentile(bv, 97.5), vc_genes=len(vc)))
T4 = pd.DataFrame(rows); T4.to_csv(os.path.join(OUT, 'E4_inheritance.csv'), index=False); K['E4'] = T4.round(4).to_dict(orient='records')
ad, ar = G.loc[grp == 'AD only', 'pct_idr'], G.loc[grp == 'AR only', 'pct_idr']
K['E4_AD_vs_AR_pct_idr_mwu_p'] = float(f'{stats.mannwhitneyu(ad, ar).pvalue:.3g}')
lm = smf.ols('pct_idr ~ AD + np.log10(length) + C(pclass)', G[grp.isin(['AD only', 'AR only'])].assign(AD=(grp == 'AD only').astype(int))).fit(cov_type='HC3')
K['E4_AD_effect_adjusted_class_length'] = [round(lm.params.AD, 3), round(lm.conf_int().loc['AD', 0], 3), round(lm.conf_int().loc['AD', 1], 3), float(f'{lm.pvalues.AD:.3g}')]

# ---- E5 is the DNA-binding-TF exception driven by the young KRAB/C2H2 zinc-finger expansion?
tf = G[(G.pclass == 'DNA-binding TFs') & G.gene_age.notna()].copy()
tf['kzfp'] = tf.family.fillna('').str.contains('Krueppel C2H2').astype(int)
e5 = {}
for lab, z in [('all DNA-binding TFs', tf), ('C2H2 zinc-finger family', tf[tf.kzfp == 1]), ('other DNA-binding TFs', tf[tf.kzfp == 0])]:
    r = stats.spearmanr(z.gene_age, z.pct_idr); med = z.groupby('age_bin').pct_idr.median().reindex(AGE); n = z.groupby('age_bin').size().reindex(AGE).fillna(0)
    m = smf.ols('gene_age ~ pct_idr', z).fit(cov_type='HC3')
    e5[lab] = dict(genes=len(z), rho=round(r[0], 3), p=float(f'{r[1]:.3g}'), median_age=float(z.gene_age.median()), median_pct_idr=round(float(z.pct_idr.median()), 2),
                   slope=[round(m.params.pct_idr, 3), round(m.conf_int().loc['pct_idr', 0], 3), round(m.conf_int().loc['pct_idr', 1], 3)],
                   median_by_age={a: (round(float(med[a]), 2) if np.isfinite(med[a]) else None) for a in AGE}, n_by_age={a: int(n[a]) for a in AGE})
m = smf.ols('gene_age ~ pct_idr * kzfp', tf).fit(cov_type='HC3')
e5['interaction_pct_idr_x_kzfp'] = [round(m.params['pct_idr:kzfp'], 3), float(f"{m.pvalues['pct_idr:kzfp']:.3g}")]
m = smf.ols('gene_age ~ pct_idr + kzfp', tf).fit(cov_type='HC3')
e5['slope_adjusted_for_kzfp'] = [round(m.params.pct_idr, 3), round(m.conf_int().loc['pct_idr', 0], 3), round(m.conf_int().loc['pct_idr', 1], 3), float(f'{m.pvalues.pct_idr:.3g}')]
K['E5_tf_exception_kzfp'] = e5

# ---- E6 age-disorder relationship within every UniProt family (even-handed; families with >= 15 genes with an age and
# genes in >= 2 age bins)
rows = []
for f, z in G.dropna(subset=['family', 'gene_age']).groupby('family'):
    if len(z) < 15 or z.age_bin.nunique() < 2: continue
    r = stats.spearmanr(z.gene_age, z.pct_idr)
    rows.append(dict(family=f, genes=len(z), pclass=z.pclass.mode().iloc[0], median_age=float(z.gene_age.median()), age_iqr=float(z.gene_age.quantile(0.75) - z.gene_age.quantile(0.25)),
                     median_pct_idr=round(float(z.pct_idr.median()), 2), rho=r[0], p=r[1]))
T6 = pd.DataFrame(rows)
if len(T6):
    from statsmodels.stats.multitest import multipletests
    T6['q'] = multipletests(T6.p.fillna(1), method='fdr_bh')[1]
    T6 = T6.sort_values('rho'); T6.to_csv(os.path.join(OUT, 'E6_age_disorder_by_family.csv'), index=False)
    K['E6_families'] = int(len(T6)); K['E6_n_positive_q05'] = int(((T6.rho > 0) & (T6.q < 0.05)).sum()); K['E6_n_negative_q05'] = int(((T6.rho < 0) & (T6.q < 0.05)).sum())
    K['E6_median_rho'] = round(float(T6.rho.median()), 3); K['E6_frac_negative'] = round(float((T6.rho < 0).mean()), 3)
    K['E6_significant'] = T6[T6.q < 0.05][['family', 'genes', 'pclass', 'rho', 'q']].round(4).to_dict(orient='records')

# between- vs within-family: slope of %IDR on gene age (per 100 Myr) without and with family fixed effects (genes in
# families with >= 15 aged genes), and with class fixed effects (all genes)
Z = G.dropna(subset=['gene_age']).copy(); Z['age100'] = Z.gene_age / 100
fams = Z.dropna(subset=['family']).groupby('family').filter(lambda z: len(z) >= 15)
def sl(m, t='age100'): return [round(m.params[t], 3), round(m.conf_int().loc[t, 0], 3), round(m.conf_int().loc[t, 1], 3), float(f'{m.pvalues[t]:.3g}')]
K['E6_slope_pct_idr_per_100Myr'] = {
    'all genes': sl(smf.ols('pct_idr ~ age100', Z).fit(cov_type='HC3')),
    'all genes, class fixed effects': sl(smf.ols('pct_idr ~ age100 + C(pclass)', Z).fit(cov_type='HC3')),
    'family genes, no family effects': sl(smf.ols('pct_idr ~ age100', fams).fit(cov_type='HC3')),
    'family genes, family fixed effects': sl(smf.ols('pct_idr ~ age100 + C(family)', fams).fit(cov_type='HC3')),
    'family genes': int(len(fams)), 'families': int(fams.family.nunique())}
fm = fams.groupby('family').agg(age=('gene_age', 'median'), idr=('pct_idr', 'median'), n=('pct_idr', 'size'))
K['E6_between_family_spearman_median_age_vs_median_idr'] = [round(v, 3) for v in stats.spearmanr(fm.age, fm.idr)] + [len(fm)]
json.dump(K, open(os.path.join(OUT, 'key_numbers.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps({k: v for k, v in K.items() if k not in ('E1', 'E2')}, indent=1, default=str))
print(T1[['pclass', 'genes', 'rho', 'p'] + [f'median_{a}' for a in AGE]].round(3).to_string()); print(T2.round(3).to_string()); print(T4.round(3).to_string())
