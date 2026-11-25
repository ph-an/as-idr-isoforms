"""Intrinsic disorder and its alternative splicing across protein classes and protein families.

Reads outputs/gene_annotations.csv (gene_annotations.py) and the isoform, exon, segment and reference-free tables.
Every analysis of the earlier TF vs non-TF comparison is repeated across ten broad protein classes (pclass), and the
disorder-relevant classifications are tested one at a time and jointly. Writes outputs/*.csv and outputs/key_numbers.json.

  A  isoform %IDR, canonical vs alternative, per class
  B  longest IDR per isoform and %IDR by gene-age bin, per class
  C  alternative vs constitutive exons: % encoding an IDR (idr_frac >= 0.5) and mean IDR fraction, per class (gene bootstrap)
  D  per-gene difference in exon IDR fraction (alternative - constitutive), per class
  E  share of a gene's unique exons that are alternative, per class
  F  exon modularity: exon boundaries inside IDR vs ordered segments relative to a same-length placement null, per class
  G  IDR density along the protein (canonical isoforms), per class
  H  IDR fates over all isoform pairs (shared genomic codons), per class
  I  UniProt protein families (>= 15 genes): disorder, isoforms, %IDR range and splicing-disorder coupling;
     variance in %IDR explained by family vs by class
  J  single classifications (localisation, function keywords, PTMs, interactors, structure, condensates, disease):
     median %IDR in vs out, and the alternative-exon disorder excess in vs out
  K  joint gene-level model of canonical %IDR on all classifications; MobiDB-lite vs MetaPredict agreement
"""
import os, json
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats

ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
HERE = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(HERE, '..'); OUT = os.path.join(HERE, 'outputs')
rng = np.random.default_rng(20261104)
G = pd.read_csv(os.path.join(OUT, 'gene_annotations.csv'), index_col=0)
CLASSES = G.groupby('pclass').pct_idr.median().sort_values().index.tolist()
K = {'classes_by_median_pct_idr': CLASSES}
AGE = ['<100', '100-500', '500-1000', '>1000']

# ---- A: isoform %IDR canonical vs alternative
iso = pd.read_csv(os.path.join(ROOT, 'data/annotated/IDRisoforms_df_geneage.csv'), low_memory=False,
                  usecols=['isoform_accession', 'base_accession', 'is_canonical', 'pct_idr', 'Length', 'max_idr_len', 'n_idr_segments', 'age_bin'])
iso['pclass'] = iso.base_accession.map(G.pclass); iso['max_idr_len'] = iso.max_idr_len.fillna(0)
rows = []
for c, d in iso.groupby('pclass'):
    can, alt = d[d.is_canonical].pct_idr, d[~d.is_canonical].pct_idr
    rows.append(dict(pclass=c, genes=d.base_accession.nunique(), canonical_n=len(can), alternative_n=len(alt), canonical_median=can.median(),
                     alternative_median=alt.median(), mwu_p=stats.mannwhitneyu(can, alt).pvalue, frac_no_idr=(d.n_idr_segments.fillna(0) == 0).mean(),
                     median_longest_idr=d[d.max_idr_len > 0].max_idr_len.median(), median_length=d.Length.median()))
TA = pd.DataFrame(rows).set_index('pclass').loc[CLASSES]; TA.to_csv(os.path.join(OUT, 'A_pct_idr_by_class.csv'))
K['A'] = TA.round(4).to_dict(orient='index')

# ---- B: age
TB = iso[iso.max_idr_len > 0].groupby(['pclass', 'age_bin']).max_idr_len.agg(['size', 'median', lambda s: s.quantile(.25), lambda s: s.quantile(.75)])
TB.columns = ['n', 'median', 'q25', 'q75']; TB.to_csv(os.path.join(OUT, 'B_longest_idr_by_class_age.csv'))
can = iso[iso.is_canonical]
K['B_spearman_age_vs_pct_idr_canonical'] = {}
ga = G.dropna(subset=['gene_age'])
for c, d in ga.groupby('pclass'):
    r = stats.spearmanr(d.gene_age, d.pct_idr); K['B_spearman_age_vs_pct_idr_canonical'][c] = [round(r[0], 3), float(f'{r[1]:.2g}'), len(d)]

# ---- C, D, E: exons
ux = pd.read_csv(os.path.join(A, 'reanalysis_2026-10-26/outputs/unique_exons.csv'))
ux['pclass'] = ux.gene.map(G.pclass); ux['age_bin'] = ux.gene.map(G.age_bin)
uidx = ux.groupby('gene').indices
def ediff(z, col='idr_frac'):
    return z[z.alt == 1][col].mean() - z[z.alt == 0][col].mean()
def boot(genes, col='idr_frac', B=500):
    pt = ediff(ux.iloc[np.concatenate([uidx[g] for g in genes])], col)
    bs = [ediff(ux.iloc[np.concatenate([uidx[g] for g in rng.choice(genes, len(genes))])], col) for _ in range(B)]
    return pt, np.percentile(bs, 2.5), np.percentile(bs, 97.5)
rows = []
for c in CLASSES:
    z = ux[ux.pclass == c]; genes = z.gene.unique()
    d, lo, hi = boot(genes); t = pd.crosstab(z.alt, z.idr_exon)
    orr = stats.contingency.odds_ratio(t.values[::-1, ::-1]).statistic if t.shape == (2, 2) else np.nan
    pg = z.groupby(['gene', 'alt']).idr_frac.mean().unstack().dropna()
    af = z.groupby('gene').alt.mean()
    rows.append(dict(pclass=c, genes=len(genes), exons_alt=int(z.alt.sum()), exons_const=int((1 - z.alt).sum()),
                     pct_idr_exon_alt=100 * z[z.alt == 1].idr_exon.mean(), pct_idr_exon_const=100 * z[z.alt == 0].idr_exon.mean(), odds_ratio=orr,
                     mean_idr_frac_alt=z[z.alt == 1].idr_frac.mean(), mean_idr_frac_const=z[z.alt == 0].idr_frac.mean(), diff=d, diff_lo=lo, diff_hi=hi,
                     paired_genes=len(pg), frac_genes_alt_higher=(pg[1] > pg[0]).mean(), frac_genes_alt_lower=(pg[1] < pg[0]).mean(),
                     wilcoxon_p=stats.wilcoxon(pg[1], pg[0]).pvalue if len(pg) > 10 else np.nan,
                     median_alt_exon_fraction=af.median(), mean_alt_exon_fraction=af.mean()))
TC = pd.DataFrame(rows).set_index('pclass'); TC.to_csv(os.path.join(OUT, 'C_exon_coupling_by_class.csv'))
K['C'] = TC.round(4).to_dict(orient='index')
pgall = ux.groupby(['gene', 'alt']).idr_frac.mean().unstack().dropna(); pgall['delta'] = pgall[1] - pgall[0]; pgall['pclass'] = pgall.index.map(G.pclass)
pgall.to_csv(os.path.join(OUT, 'D_gene_delta.csv'))
afr = ux.groupby('gene').alt.mean().rename('alt_fraction').to_frame(); afr['pclass'] = afr.index.map(G.pclass)
afr.to_csv(os.path.join(OUT, 'E_gene_alt_fraction.csv'))
K['E_kruskal_alt_fraction_across_classes_p'] = float(f'{stats.kruskal(*[g.alt_fraction for _, g in afr.groupby("pclass")]).pvalue:.3g}')
# coupling by class x age
rows = []
for c in CLASSES:
    for a in AGE:
        genes = ux[(ux.pclass == c) & (ux.age_bin == a)].gene.unique()
        if len(genes) < 30: continue
        d, lo, hi = boot(genes, B=300)
        rows.append(dict(pclass=c, age_bin=a, genes=len(genes), diff=d, lo=lo, hi=hi,
                         pct_alt_idr=100 * ux[(ux.gene.isin(genes)) & (ux.alt == 1)].idr_exon.mean()))
TCA = pd.DataFrame(rows); TCA.to_csv(os.path.join(OUT, 'C_coupling_by_class_age.csv'))
rows = []
for a in AGE:
    genes = ux[ux.age_bin == a].gene.unique(); d, lo, hi = boot(genes, B=500)
    rows.append(dict(age_bin=a, genes=len(genes), diff=d, lo=lo, hi=hi))
K['C_coupling_by_age_all'] = pd.DataFrame(rows).set_index('age_bin').round(4).to_dict(orient='index')

# ---- F: modularity by class
S = pd.read_csv(os.path.join(A, 'reanalysis_2026-10-26/outputs/R6_segments_obs_exp.csv'))
S['gene'] = S.iso.str.split('-').str[0]; S['pclass'] = S.gene.map(G.pclass)
TF_ = S.groupby(['pclass', 'kind']).agg(n=('obs', 'size'), obs=('obs', 'sum'), exp=('exp', 'sum'), single=('obs', lambda x: (x == 0).mean()))
TF_['O_E'] = TF_.obs / TF_.exp; TF_.to_csv(os.path.join(OUT, 'F_modularity_by_class.csv'))
K['F_OE'] = TF_.O_E.unstack().round(3).to_dict(orient='index')

# ---- G: positional density (canonical isoforms)
idr = pd.read_csv(os.path.join(ROOT, 'data/annotated/idr_df.csv')); idr = idr[idr.idr_aaLen >= 20]
cs = can.set_index('isoform_accession'); prof = {}
for c in CLASSES:
    accs = cs[cs.pclass == c]; acc = np.zeros(50); n = 0
    sel = idr[idr.isoform_accession.isin(accs.index)]
    lens = accs.Length.to_dict()
    for a, d in sel.groupby('isoform_accession'):
        L = lens[a]; m = np.zeros(L, bool)
        for s, e in zip(d.idr_start_aa, d.idr_end_aa): m[int(s) - 1:int(e)] = True
        acc += np.array([b.mean() for b in np.array_split(m, 50)])
    n = len(accs); prof[c] = (acc / n).round(4).tolist()
json.dump(prof, open(os.path.join(OUT, 'G_position_profiles.json'), 'w'))
K['G_terminal_vs_middle'] = {c: dict(ends=round(np.mean(p[:5] + p[-5:]), 3), middle=round(np.mean(p[20:30]), 3)) for c, p in prof.items()}

# ---- H: fates (codon-level, all pairs)
F = pd.read_csv(os.path.join(A, 'reference_free/outputs/all_pairs_idr_fates.csv')); F['pclass'] = F.gene.map(G.pclass)
FT = ['intact', 'internal deletion', 'truncated', 'removed']
TH = F.pivot_table(index='pclass', columns='fate', values='n', aggfunc='sum').reindex(index=CLASSES, columns=FT).fillna(0)
TH = 100 * TH.div(TH.sum(axis=1), axis=0); TH.to_csv(os.path.join(OUT, 'H_fates_by_class.csv')); K['H'] = TH.round(1).to_dict(orient='index')
RF = pd.read_csv(os.path.join(A, 'reference_free/outputs/gene_core_variable.csv')).set_index('gene'); RF['pclass'] = RF.index.map(G.pclass)
K['H_core_vs_variable_by_class'] = RF.groupby('pclass')[['dis_core', 'dis_var']].mean().round(3).to_dict(orient='index')

# ---- I: families
fam = G.dropna(subset=['family']).copy()
big = fam.family.value_counts(); big = big[big >= 15].index
rows = []
for f in big:
    d = fam[fam.family == f]; z = ux[ux.gene.isin(d.index)]
    rows.append(dict(family=f, genes=len(d), pclass=d.pclass.mode().iloc[0], median_pct_idr=d.pct_idr.median(), q25=d.pct_idr.quantile(.25), q75=d.pct_idr.quantile(.75),
                     mean_iso=d.n_iso.mean(), median_range=d.idr_range.median(), frac_tf=d.is_tf.mean(), median_age=d.gene_age.median(),
                     coupling=ediff(z) if z.alt.nunique() == 2 else np.nan, exon_genes=z.gene.nunique(), mobidb_median=d.mobidb_pct.median()))
TI = pd.DataFrame(rows).sort_values('median_pct_idr', ascending=False); TI.to_csv(os.path.join(OUT, 'I_families.csv'), index=False)
K['I_n_families_ge15'] = len(TI); K['I_genes_in_them'] = int(TI.genes.sum())
fb = fam[fam.family.isin(big)]
def eta2(y, g):
    m = y.mean(); return ((y.groupby(g).transform('mean') - m) ** 2).sum() / ((y - m) ** 2).sum()
e_f = eta2(fb.pct_idr, fb.family); e_c = eta2(fb.pct_idr, fb.pclass)
null = [eta2(fb.pct_idr, pd.Series(rng.permutation(fb.family.values), index=fb.index)) for _ in range(200)]
K['I_variance_explained_pct_idr'] = dict(genes=len(fb), family=round(e_f, 3), pclass=round(e_c, 3), family_null_mean=round(np.mean(null), 3),
                                         family_null_95=round(np.percentile(null, 95), 3))
e_fr = eta2(fb.idr_range, fb.family); null_r = [eta2(fb.idr_range, pd.Series(rng.permutation(fb.family.values), index=fb.index)) for _ in range(200)]
K['I_variance_explained_idr_range'] = dict(family=round(e_fr, 3), family_null_mean=round(np.mean(null_r), 3), family_null_95=round(np.percentile(null_r, 95), 3))
K['I_top10'] = TI.head(10)[['family', 'genes', 'median_pct_idr', 'coupling']].round(3).to_dict(orient='records')
K['I_bottom10'] = TI.tail(10)[['family', 'genes', 'median_pct_idr', 'coupling']].round(3).to_dict(orient='records')
tc = TI.dropna(subset=['coupling'])
K['I_family_coupling_vs_disorder_spearman'] = [round(x, 3) for x in stats.spearmanr(tc.median_pct_idr, tc.coupling)] + [len(tc)]

# ---- J: single classifications
BIN = [c for c in G.columns if c.startswith(('loc_', 'kw_'))] + ['enzyme', 'transmembrane', 'signal_peptide', 'disease', 'has_pdb', 'disprot',
                                                                 'condensate', 'coiled_coil']
G['many_phospho'] = (G.n_phospho >= 5).astype(float); G['many_interactors'] = (G.n_interactors >= 10).astype(float)
BIN += ['many_phospho', 'many_interactors']
ux_g = ux.set_index('gene')
rows = []
for b in BIN:
    x = G[b].dropna().astype(bool); inn, out = x.index[x], x.index[~x]
    if len(inn) < 50: continue
    zi, zo = ux[ux.gene.isin(inn)], ux[ux.gene.isin(out)]
    rows.append(dict(feature=b, genes_in=len(inn), median_pct_idr_in=G.loc[inn, 'pct_idr'].median(), median_pct_idr_out=G.loc[out, 'pct_idr'].median(),
                     mwu_p=stats.mannwhitneyu(G.loc[inn, 'pct_idr'], G.loc[out, 'pct_idr']).pvalue,
                     mean_iso_in=G.loc[inn, 'n_iso'].mean(), mean_iso_out=G.loc[out, 'n_iso'].mean(),
                     coupling_in=ediff(zi), coupling_out=ediff(zo)))
TJ = pd.DataFrame(rows)
# exon-level interaction test of coupling difference (gene-clustered)
pv = []
for r in TJ.itertuples():
    z = ux[['gene', 'idr_frac', 'alt', 'log_len', 'pos', 'log_ntx']].copy(); z['f'] = z.gene.map(G[r.feature].astype(float)).fillna(0)
    m = smf.ols('idr_frac ~ alt * f + log_len + C(pos) + log_ntx', z).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(z.gene)[0]})
    pv.append((m.params['alt:f'], m.pvalues['alt:f']))
TJ['coupling_interaction'] = [p[0] for p in pv]; TJ['coupling_interaction_p'] = [p[1] for p in pv]
from statsmodels.stats.multitest import multipletests
TJ['mwu_q'] = multipletests(TJ.mwu_p, method='fdr_bh')[1]; TJ['coupling_q'] = multipletests(TJ.coupling_interaction_p, method='fdr_bh')[1]
TJ = TJ.sort_values('median_pct_idr_in', ascending=False); TJ.to_csv(os.path.join(OUT, 'J_features.csv'), index=False)
K['J'] = TJ.round(4).to_dict(orient='records')

# ---- K: joint model and MobiDB agreement
J = G.copy(); J['log_len'] = np.log10(J.length)
for b in BIN: J[b] = J[b].astype(float)
J = J.dropna(subset=[b for b in BIN if b.startswith('loc_')])
form = 'pct_idr ~ log_len + ' + ' + '.join(BIN)
mk = smf.ols(form, J).fit(cov_type='HC3')
TK = pd.DataFrame({'coef': mk.params, 'lo': mk.conf_int()[0], 'hi': mk.conf_int()[1], 'p': mk.pvalues}).drop('Intercept')
TK.to_csv(os.path.join(OUT, 'K_joint_model.csv')); K['K_joint_model'] = dict(n=int(mk.nobs), r2=round(mk.rsquared, 3))
fm = smf.ols('pct_idr ~ log_len + C(pclass)', J).fit(); K['K_r2_class_only'] = round(fm.rsquared, 3)
mb = G.dropna(subset=['mobidb_pct'])
K['K_mobidb_vs_metapredict'] = dict(genes=len(mb), spearman=round(stats.spearmanr(mb.mobidb_pct, mb.pct_idr)[0], 3),
                                    class_medians={c: [round(d.pct_idr.median(), 1), round(d.mobidb_pct.median(), 1)] for c, d in mb.groupby('pclass')})
json.dump(K, open(os.path.join(OUT, 'key_numbers.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps({k: K[k] for k in ['classes_by_median_pct_idr', 'I_variance_explained_pct_idr', 'I_variance_explained_idr_range', 'I_family_coupling_vs_disorder_spearman', 'K_joint_model', 'K_r2_class_only', 'C_coupling_by_age_all']}, indent=1, default=str))
print(TC[['genes', 'pct_idr_exon_alt', 'pct_idr_exon_const', 'diff', 'diff_lo', 'diff_hi', 'frac_genes_alt_higher', 'median_alt_exon_fraction']].round(3).to_string())
print(TH.round(1).to_string())
print(TJ[['feature', 'genes_in', 'median_pct_idr_in', 'median_pct_idr_out', 'coupling_in', 'coupling_out', 'coupling_q']].round(3).to_string())
