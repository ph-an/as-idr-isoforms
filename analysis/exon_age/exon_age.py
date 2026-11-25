"""Exon age from splice-site conservation (Minkin & Salzberg 2025, NAR 53:gkaf184) and its relation to disorder and splicing.

Input: data/raw/splice_sites_conservation/splice_sites.csv.gz (one row per splice site per transcript; cons_GTAG = number of
the 405 aligned mammalian genomes in which the GT donor or AG acceptor dinucleotide is conserved). Decoy ('Random') sites are
dropped; a site's conservation is its maximum over annotations.

Exon conservation: each coding exon (gene, Ensembl exon) is matched by GRCh38 coordinate to its acceptor (start - 2 on +,
end + 2 on -) and donor (end + 1 on +, start - 1 on -). Internal exons with both sites take the lower count, because an exon
is no older than its youngest boundary. Bins: lineage-restricted (< 100 of 405 genomes), intermediate (100-349) and deeply
conserved (>= 350). The count is a conservation score, not a date: it depends on how densely each clade is sampled.

X1  exon disorder by conservation, all exons and by usage
X2  within-gene models (gene fixed effects; gene-clustered errors): exon IDR fraction ~ conservation + log length + usage + position
X3  alternative usage and inclusion by conservation, within genes
X4  composition of disordered residues from lineage-restricted vs deeply conserved exons of the same gene (gene-paired)
X5  by protein class and gene age (GenOrigin): does the young-exon effect hold inside old genes?
"""
import os, json
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); A = os.path.join(ROOT, 'analysis')
OUT = os.path.join(A, 'exon_age', 'outputs'); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261104)
K = {}

# ---- sites
sp = os.path.join(OUT, 'sites_unique.parquet')
if not os.path.exists(sp):
    S = pd.read_csv(os.path.join(ROOT, 'data/raw/splice_sites_conservation/splice_sites.csv.gz'),
                    usecols=['dataset', 'site_type', 'chr', 'strand', 'pos', 'cons_GTAG', 'motif', 'phastCons_0', 'phastCons_1'])
    S = S[S.dataset != 'Random']
    S.groupby(['chr', 'strand', 'pos', 'site_type']).agg(cons=('cons_GTAG', 'max'), motif=('motif', 'first'),
                                                         phast=('phastCons_0', 'mean')).reset_index().to_parquet(sp)
U = pd.read_parquet(sp); key = U.set_index(['chr', 'strand', 'pos', 'site_type']).cons
# alignability: median conservation of the 30 exonic bases next to the site (acceptor: downstream; donor: upstream), from the
# per-position columns; written into sites_unique.parquet by the chunked pass below if missing
if 'flank_neg' not in U.columns:
    neg = [f'cons_{i}' for i in range(-30, 0)]; posc = [f'cons_{i}' for i in range(2, 32)]; parts = []
    for ch in pd.read_csv(os.path.join(ROOT, 'data/raw/splice_sites_conservation/splice_sites.csv.gz'), usecols=['dataset', 'site_type', 'chr', 'strand', 'pos'] + neg + posc,
                          chunksize=200000, dtype={c: 'float32' for c in neg + posc}):
        ch = ch[ch.dataset != 'Random']; o = ch[['site_type', 'chr', 'strand', 'pos']].copy()
        o['flank_neg'] = ch[neg].median(axis=1); o['flank_pos'] = ch[posc].median(axis=1); parts.append(o)
    F = pd.concat(parts).groupby(['chr', 'strand', 'pos', 'site_type'])[['flank_neg', 'flank_pos']].max().reset_index()
    U = U.merge(F, on=['chr', 'strand', 'pos', 'site_type'], how='left'); U.to_parquet(sp)
U['exon_flank'] = np.where(U.site_type == 'a', U.flank_pos, U.flank_neg)
fkey = U.set_index(['chr', 'strand', 'pos', 'site_type']).exon_flank
K['sites'] = len(U)

# ---- exons
E = pd.read_parquet(os.path.join(ROOT, 'data/genomic/processed/exon_level_table.parquet'),
                    columns=['isoform_accession', 'base_accession', 'exon_id', 'chr', 'start', 'end', 'strand', 'is_coding_exon', 'aa_start', 'aa_end',
                             'transcript_inclusion_fraction'])
E = E[E.is_coding_exon]
X = E.drop_duplicates(['base_accession', 'exon_id']).copy()
p = (X.strand == '+').values
acc = np.where(p, X.start - 2, X.end + 2); don = np.where(p, X.end + 1, X.start - 1)
look = lambda pos, t: key.reindex(pd.MultiIndex.from_arrays([X.chr.values, X.strand.values, pos, np.array([t] * len(X))])).values
X['cons_acc'] = look(acc, 'a'); X['cons_don'] = look(don, 'd')
flook = lambda pos, t: fkey.reindex(pd.MultiIndex.from_arrays([X.chr.values, X.strand.values, pos, np.array([t] * len(X))])).values
X['flank'] = np.fmin(flook(acc, 'a'), flook(don, 'd'))
X['n_sites'] = np.isfinite(X.cons_acc).astype(int) + np.isfinite(X.cons_don).astype(int)
X['cons'] = X[['cons_acc', 'cons_don']].min(axis=1)
incl = E.groupby(['base_accession', 'exon_id']).transcript_inclusion_fraction.mean()
X['inclusion'] = incl.reindex(pd.MultiIndex.from_arrays([X.base_accession, X.exon_id])).values
UX = pd.read_csv(os.path.join(A, 'reanalysis_2026-10-26/outputs/unique_exons.csv'))
D = UX.merge(X[['base_accession', 'exon_id', 'cons_acc', 'cons_don', 'n_sites', 'cons', 'inclusion', 'flank']], left_on=['gene', 'exon_id'], right_on=['base_accession', 'exon_id'], how='left')
G = pd.read_csv(os.path.join(A, 'protein_classes/outputs/gene_annotations.csv'), index_col=0)
D['pclass'] = D.gene.map(G.pclass); D['gene_age'] = D.gene.map(G.gene_age); D['age_bin'] = D.gene.map(G.age_bin)
K['unique_exons'] = len(D); K['with_any_site'] = int((D.n_sites > 0).sum()); K['internal_with_both_sites'] = int(((D.pos == 'internal') & (D.n_sites == 2)).sum())
BINS = [(-1, 99, 'lineage-restricted (<100)'), (99, 349, 'intermediate (100-349)'), (349, 405, 'deeply conserved (>=350)')]
LAB = [b[2] for b in BINS]
D['cbin'] = pd.cut(D.cons, [-1, 99, 349, 405], labels=LAB)
D['young'] = (D.cons < 100).astype(float); D.loc[D.cons.isna(), 'young'] = np.nan
D['cons100'] = D.cons / 100
D.to_csv(os.path.join(OUT, 'exon_conservation.csv'), index=False)
I = D[(D.pos == 'internal') & (D.n_sites == 2)].copy()      # main set: internal exons with both sites
K['cons_quantiles_internal'] = I.cons.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).round(0).to_dict()
K['bin_sizes_internal'] = I.cbin.value_counts().reindex(LAB).astype(int).to_dict()

# ---- X1
def summ(z):
    return dict(exons=len(z), genes=z.gene.nunique(), mean_idr=round(z.idr_frac.mean(), 4), pct_idr_exon=round(100 * (z.idr_frac >= 0.5).mean(), 2),
                pct_alt=round(100 * z.alt.mean(), 2), median_len_aa=float(z.aa_len.median()))
X1 = pd.DataFrame({b: summ(I[I.cbin == b]) for b in LAB}).T
X1u = pd.DataFrame({(b, u): summ(I[(I.cbin == b) & (I.alt == u)]) for b in LAB for u in (0, 1)}).T
X1.to_csv(os.path.join(OUT, 'X1_disorder_by_conservation.csv')); X1u.to_csv(os.path.join(OUT, 'X1_disorder_by_conservation_usage.csv'))
K['X1'] = X1.to_dict(orient='index'); K['X1_by_usage'] = {f'{b}|{"alt" if u else "const"}': v for (b, u), v in X1u.to_dict(orient='index').items()}
# first / last exons (one site) reported separately
F = D[D.pos.isin(['first', 'last']) & (D.n_sites >= 1)]
K['X1_terminal_exons'] = {f'{ps}|{b}': summ(F[(F.pos == ps) & (F.cbin == b)]) for ps in ('first', 'last') for b in LAB}

# ---- X2 within-gene
def fe(df, y, xs):
    z = df.dropna(subset=[y] + xs).copy(); z = z[z.groupby('gene')[y].transform('size') > 1]
    for c in [y] + xs: z[c + '_d'] = z[c] - z.groupby('gene')[c].transform('mean')
    m = smf.ols(f'{y}_d ~ ' + ' + '.join(c + '_d' for c in xs) + ' - 1', z).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(z.gene)[0]})
    return m, z
I['first_last'] = 0
m, z = fe(I, 'idr_frac', ['young', 'log_len', 'alt', 'log_ntx'])
K['X2_fe_young'] = dict(genes=int(z.gene.nunique()), exons=len(z), coef=round(m.params.young_d, 4), lo=round(m.conf_int().loc['young_d', 0], 4),
                        hi=round(m.conf_int().loc['young_d', 1], 4), p=float(f'{m.pvalues.young_d:.3g}'), alt_coef=round(m.params.alt_d, 4))
m2, z2 = fe(I, 'idr_frac', ['cons100', 'log_len', 'alt', 'log_ntx'])
K['X2_fe_per100genomes'] = [round(m2.params.cons100_d, 4), round(m2.conf_int().loc['cons100_d', 0], 4), round(m2.conf_int().loc['cons100_d', 1], 4), float(f'{m2.pvalues.cons100_d:.3g}')]
for u, lab in [(0, 'constitutive'), (1, 'alternative')]:
    mu, zu = fe(I[I.alt == u], 'idr_frac', ['young', 'log_len', 'log_ntx'])
    K[f'X2_fe_young_{lab}_only'] = [round(mu.params.young_d, 4), round(mu.conf_int().loc['young_d', 0], 4), round(mu.conf_int().loc['young_d', 1], 4), float(f'{mu.pvalues.young_d:.3g}'), int(zu.gene.nunique())]
# alignment control: low GT/AG counts could reflect exons that fail to align rather than young splice sites. (a) add the
# exonic-flank conservation as a covariate; (b) keep only exons whose flanking exonic bases align in >= 200 genomes, so that a
# low count means the splice signal is absent where the exon itself aligns
I['flank100'] = I.flank / 100
ma, za = fe(I, 'idr_frac', ['young', 'flank100', 'log_len', 'alt', 'log_ntx'])
K['X2_fe_young_adj_flank'] = [round(ma.params.young_d, 4), round(ma.conf_int().loc['young_d', 0], 4), round(ma.conf_int().loc['young_d', 1], 4), float(f'{ma.pvalues.young_d:.3g}')]
Ia = I[I.flank >= 200]
mb_, zb_ = fe(Ia, 'idr_frac', ['young', 'log_len', 'alt', 'log_ntx'])
K['X2_fe_young_alignable_only'] = dict(young_exons=int(Ia.young.sum()), genes=int(zb_.gene.nunique()), coef=round(mb_.params.young_d, 4),
                                       lo=round(mb_.conf_int().loc['young_d', 0], 4), hi=round(mb_.conf_int().loc['young_d', 1], 4), p=float(f'{mb_.pvalues.young_d:.3g}'),
                                       mean_idr_young=round(Ia[Ia.young == 1].idr_frac.mean(), 4), mean_idr_old=round(Ia[Ia.cons >= 350].idr_frac.mean(), 4))
rob = {}
for t in (200, 250, 300, 350):
    mt_, zt_ = fe(I[I.flank >= t], 'idr_frac', ['young', 'log_len', 'alt', 'log_ntx'])
    rob[t] = [round(mt_.params.young_d, 4), round(mt_.conf_int().loc['young_d', 0], 4), round(mt_.conf_int().loc['young_d', 1], 4), float(f'{mt_.pvalues.young_d:.3g}'), int(I[I.flank >= t].young.sum())]
K['X2_fe_young_by_flank_threshold'] = rob
I['fbin'] = pd.qcut(I.flank, 10, duplicates='drop')
K['X2_flank_stratified_mean_idr'] = {str(k): v for k, v in I.groupby(['fbin', 'young'], observed=True).idr_frac.mean().unstack().round(4).to_dict(orient='index').items()}
K['X2_flank_by_bin'] = I.groupby('cbin', observed=True).flank.median().to_dict()
# paired genes: genes with both lineage-restricted and deeply conserved internal exons
pg = I[I.cbin.isin([LAB[0], LAB[2]])].groupby(['gene', 'cbin'], observed=True).idr_frac.mean().unstack()
pg = pg.dropna()
K['X2_paired_genes'] = dict(genes=len(pg), mean_young=round(pg[LAB[0]].mean(), 4), mean_old=round(pg[LAB[2]].mean(), 4),
                            frac_young_higher=round((pg[LAB[0]] > pg[LAB[2]]).mean(), 3), frac_young_lower=round((pg[LAB[0]] < pg[LAB[2]]).mean(), 3),
                            wilcoxon_p=float(f'{stats.wilcoxon(pg[LAB[0]], pg[LAB[2]]).pvalue:.3g}'))
# dose: idr by finer conservation deciles (internal exons)
I['cdec'] = pd.cut(I.cons, [-1, 19, 49, 99, 199, 299, 349, 379, 395, 405])
X2d = I.groupby('cdec', observed=True).agg(exons=('idr_frac', 'size'), mean_idr=('idr_frac', 'mean'), pct_alt=('alt', 'mean'), inclusion=('inclusion', 'mean'))
X2d.to_csv(os.path.join(OUT, 'X2_by_conservation_decile.csv')); K['X2_deciles'] = {str(k): v for k, v in X2d.round(4).to_dict(orient='index').items()}

# ---- X3 usage
m3, z3 = fe(I, 'alt', ['young', 'log_len', 'log_ntx'])
K['X3_fe_alt_on_young'] = [round(m3.params.young_d, 4), round(m3.conf_int().loc['young_d', 0], 4), round(m3.conf_int().loc['young_d', 1], 4), float(f'{m3.pvalues.young_d:.3g}')]
K['X3_inclusion_by_bin'] = I.groupby('cbin', observed=True).inclusion.mean().round(4).to_dict()
K['X3_spearman_cons_inclusion_alt_exons'] = [round(v, 4) for v in stats.spearmanr(I[I.alt == 1].cons, I[I.alt == 1].inclusion, nan_policy='omit')]

# ---- X4 composition of disordered residues by exon conservation (gene-paired)
iso = pd.read_csv(os.path.join(ROOT, 'data/annotated/IDRisoforms_df.csv'), usecols=['isoform_accession', 'Sequence']); seq = dict(zip(iso.isoform_accession, iso.Sequence))
idr = pd.read_csv(os.path.join(ROOT, 'data/annotated/idr_df.csv')); idr = idr[idr.idr_aaLen >= 20]
segs = {a: [(int(s), int(e)) for s, e in zip(d.idr_start_aa, d.idr_end_aa)] for a, d in idr.groupby('isoform_accession')}
occ = E.dropna(subset=['aa_start']).drop_duplicates(['base_accession', 'exon_id']).set_index(['base_accession', 'exon_id'])
AA = {'aro': 'FWY', 'acid': 'DE', 'basic': 'KR', 'his': 'H', 'cys': 'C', 'pro': 'P', 'st': 'ST', 'qn': 'QN', 'gly': 'G', 'hyd': 'AILMV'}
pool = {}
for r in I[I.cbin.isin([LAB[0], LAB[2]])].itertuples():
    try: o = occ.loc[(r.gene, r.exon_id)]
    except KeyError: continue
    a = o.isoform_accession
    if a not in seq: continue
    s = seq[a]; lo, hi = int(o.aa_start), int(o.aa_end)
    m = np.zeros(len(s) + 1, bool)
    for st, en in segs.get(a, []): m[st:en + 1] = True
    res = ''.join(s[i - 1] for i in range(lo, min(hi, len(s)) + 1) if m[i])
    if res: pool[(r.gene, r.cbin)] = pool.get((r.gene, r.cbin), '') + res
rows = []
for (g, b), s in pool.items():
    if len(s) < 10: continue
    d = {k: sum(s.count(c) for c in v) / len(s) for k, v in AA.items()}
    d['ncpr'] = d['basic'] - d['acid']; d['fcr'] = d['basic'] + d['acid']; d.update(gene=g, cbin=b, n=len(s)); rows.append(d)
C = pd.DataFrame(rows)
W = C.pivot(index='gene', columns='cbin', values=list(AA) + ['ncpr', 'fcr']).dropna()
x4 = {}
for f in list(AA) + ['ncpr', 'fcr']:
    d = 100 * (W[(f, LAB[0])] - W[(f, LAB[2])]).values
    bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(1000)]
    x4[f] = dict(diff_pp=round(d.mean(), 3), lo=round(np.percentile(bs, 2.5), 3), hi=round(np.percentile(bs, 97.5), 3), p=float(f'{stats.wilcoxon(d).pvalue:.3g}'))
from statsmodels.stats.multitest import multipletests
q = multipletests([v['p'] for v in x4.values()], method='fdr_bh')[1]
for (f, v), qq in zip(x4.items(), q): v['q'] = float(f'{qq:.3g}')
K['X4_genes'] = len(W); K['X4_young_minus_old_disorder_composition'] = x4
pd.DataFrame(x4).T.to_csv(os.path.join(OUT, 'X4_composition_young_vs_old_exons.csv'))
# compare with gene-age and terminal-minority signatures
ev = json.load(open(os.path.join(A, 'evolution/outputs/key_numbers.json')))['E3_young_minus_old_core_pp']
rf = pd.read_csv(os.path.join(A, 'reference_free/outputs/composition_vs_core.csv'))
mt = rf[rf.kind == 'minority|terminal'].set_index('feature').diff_pp; vi = rf[rf.kind == 'variable|internal'].set_index('feature').diff_pp
common = [f for f in x4 if f in ev and f in mt.index and f != 'fcr']
xv = [x4[f]['diff_pp'] for f in common]
K['X4_similarity_spearman'] = {'vs_young_gene_core': [round(v, 3) for v in stats.spearmanr(xv, [ev[f] for f in common])],
                               'vs_minority_terminal_variable': [round(v, 3) for v in stats.spearmanr(xv, [mt[f] for f in common])],
                               'vs_internal_variable': [round(v, 3) for v in stats.spearmanr(xv, [vi[f] for f in common])], 'features': common}

# ---- X5 class and gene age
rows = []
for c, z in I.groupby('pclass'):
    try:
        mc, zc = fe(z, 'idr_frac', ['young', 'log_len', 'alt', 'log_ntx'])
        rows.append(dict(pclass=c, genes=zc.gene.nunique(), pct_young=round(100 * z.young.mean(), 2), coef=mc.params.young_d,
                         lo=mc.conf_int().loc['young_d', 0], hi=mc.conf_int().loc['young_d', 1], p=mc.pvalues.young_d))
    except Exception: pass
X5 = pd.DataFrame(rows); X5.to_csv(os.path.join(OUT, 'X5_young_exon_effect_by_class.csv'), index=False); K['X5_by_class'] = X5.round(4).to_dict(orient='records')
rows = []
for b in ['<100', '100-500', '500-1000', '>1000']:
    z = I[I.age_bin == b]
    mb, zb = fe(z, 'idr_frac', ['young', 'log_len', 'alt', 'log_ntx'])
    rows.append(dict(age_bin=b, genes=zb.gene.nunique(), pct_young_exons=round(100 * z.young.mean(), 2), median_cons=float(z.cons.median()),
                     coef=mb.params.young_d, lo=mb.conf_int().loc['young_d', 0], hi=mb.conf_int().loc['young_d', 1], p=mb.pvalues.young_d))
X5a = pd.DataFrame(rows); X5a.to_csv(os.path.join(OUT, 'X5_young_exon_effect_by_gene_age.csv'), index=False); K['X5_by_gene_age'] = X5a.round(4).to_dict(orient='records')
gm = I.groupby('gene').agg(med_cons=('cons', 'median'), age=('gene_age', 'first')).dropna()
K['X5_gene_median_cons_vs_gene_age_spearman'] = [round(v, 4) for v in stats.spearmanr(gm.med_cons, gm.age)]


# ---- X6 exon IDR fraction by splice-site conservation within each UniProt family (families with >= 15 genes and >= 10
# internal exons in every conservation bin); within-family difference = lineage-restricted/intermediate (<350) minus deeply
# conserved (>=380) exons, with a gene bootstrap
I['family'] = I.gene.map(G.family)
FB = [(-1, 349, '<350'), (349, 379, '350-379'), (379, 395, '380-395'), (395, 405, '396-405')]
I['fbin6'] = pd.cut(I.cons, [b[0] for b in FB] + [405], labels=[b[2] for b in FB])
big = G.family.value_counts(); big = big[big >= 15].index
rows = []
for f, z in I[I.family.isin(big)].groupby('family'):
    n = z.fbin6.value_counts()
    if n.min() < 10: continue
    m = z.groupby('fbin6', observed=False).idr_frac.mean()
    young = z[z.cons < 350]; old = z[z.cons >= 380]
    gidx = {g: d.index for g, d in z.groupby('gene')}; genes = np.array(list(gidx))
    bs = []
    for _ in range(500):
        zz = z.loc[np.concatenate([gidx[g] for g in rng.choice(genes, len(genes))])]
        a, b = zz[zz.cons < 350].idr_frac, zz[zz.cons >= 380].idr_frac
        if len(a) and len(b): bs.append(a.mean() - b.mean())
    rows.append(dict(family=f, pclass=G[G.family == f].pclass.mode().iloc[0], genes=z.gene.nunique(), exons=len(z), **{f'mean_{k}': m[k] for k in m.index},
                     **{f'n_{k}': int(n[k]) for k in n.index}, diff=young.idr_frac.mean() - old.idr_frac.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5)))
X6 = pd.DataFrame(rows).sort_values('diff'); X6.to_csv(os.path.join(OUT, 'X6_exon_idr_by_conservation_by_family.csv'), index=False)
K['X6_families'] = len(X6); K['X6_n_positive'] = int((X6['diff'] > 0).sum()); K['X6_n_ci_above_0'] = int((X6.lo > 0).sum()); K['X6_n_ci_below_0'] = int((X6.hi < 0).sum())
K['X6_median_diff'] = round(float(X6['diff'].median()), 4)
json.dump(K, open(os.path.join(OUT, 'key_numbers.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps({k: v for k, v in K.items() if k not in ('X2_deciles',)}, indent=1, default=str)[:9000])
print(X2d.round(3).to_string())
