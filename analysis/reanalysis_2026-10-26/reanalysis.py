"""Exon usage and predicted disorder: unique-exon table, sensitivity to the exon-usage definition, exon boundaries.

Reads the committed pipeline outputs (no raw inputs are re-processed) and writes outputs/unique_exons.csv (read by the
other analyses), outputs/key_numbers.json and the exon-boundary tables.

Inputs (paths relative to the repository root; override with TFIDR_ROOT):
  data/genomic/analysis/exon_idr_annotation.csv     s9 per (isoform x coding exon) table with mapped-transcript usage
  data/genomic/processed/exon_level_table.parquet   exon table (TF family, GENCODE-wide inclusion)
  data/genomic/interim/{idr,nonidr}_exon_overlap_table.parquet

Analyses
  R5  sensitivity: GENCODE-wide alternative definition and continuous inclusion fraction.
  R6  exon-boundary depletion within IDRs against a same-protein position-randomisation null (length-controlled).
Runtime ~5 min on a laptop.
"""
import os, json
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats

ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
GEN = os.path.join(ROOT, 'data', 'genomic')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')
os.makedirs(OUT, exist_ok=True)
K = {}  # key numbers
cl = lambda d: dict(cov_type='cluster', cov_kwds={'groups': pd.factorize(d.gene)[0]})   # gene-clustered errors

# ---------------------------------------------------------------- load
A = pd.read_csv(os.path.join(GEN, 'analysis', 'exon_idr_annotation.csv'))
E = pd.read_parquet(os.path.join(GEN, 'processed', 'exon_level_table.parquet'))
fam = E.drop_duplicates('base_accession').set_index('base_accession').tf_family
A['fam'] = A.gene.map(fam).astype(str)
A['grp'] = np.where(~A.is_tf, 'nonTF', np.where(A.fam == 'C2H2 ZF', 'C2H2', 'TF_other'))
AU = A[A.usage.isin(['alternative', 'constitutive'])].copy()
AU['aa_len'] = AU.aa_end - AU.aa_start + 1
AU['first'] = AU.exon_rank == 1
AU['last'] = AU.exon_rank == AU.n_coding_exons
g = (AU.groupby(['gene', 'exon_id'])
       .agg(idr_frac=('idr_frac', 'mean'), usage=('usage', 'first'), is_tf=('is_tf', 'first'), grp=('grp', 'first'),
            aa_len=('aa_len', 'median'), first=('first', 'max'), last=('last', 'max'), n_tx=('n_tx_gene', 'first'))
       .reset_index())
g['alt'] = (g.usage == 'alternative').astype(int); g['tf'] = g.is_tf.astype(int)
g['log_len'] = np.log10(g.aa_len.clip(lower=1)); g['log_ntx'] = np.log2(g.n_tx)
g['pos'] = np.select([g['first'], g['last']], ['first', 'last'], 'internal')
g['idr_exon'] = (g.idr_frac >= 0.5).astype(int)
g = g.join(g.groupby('gene').apply(lambda d: np.average(d.idr_frac, weights=d.aa_len), include_groups=False)
            .rename('gene_idr'), on='gene')
g.to_csv(os.path.join(OUT, 'unique_exons.csv'), index=False)
K['n_unique_exons'] = int(len(g)); K['n_genes_ge2_tx'] = int(g.gene.nunique())


# ---------------------------------------------------------------- R5 sensitivity to the alternative-exon definition
Ec = E[E.is_coding_exon & E.idr_frac.notna()]
u = Ec.groupby(['gene_id', 'exon_id']).agg(idr=('idr_frac', 'mean'), tf=('is_tf', 'first'), gw_const=('is_constitutive', 'first'),
                                           incl=('transcript_inclusion_fraction', 'first'), ntx=('n_coding_transcripts_gene', 'first')).reset_index()
u = u[u.ntx >= 2].rename(columns={'gene_id': 'gene'}); u['alt'] = (~u.gw_const).astype(int); u['tfi'] = u.tf.astype(int)
m5 = smf.ols('idr ~ alt*tfi', u).fit(**cl(u))
m5c = smf.ols('idr ~ incl*tfi + np.log2(ntx)', u).fit(**cl(u))
K['R5_gencode_wide'] = {t: [round(m5.params[t], 4), round(m5.conf_int().loc[t, 0], 4), round(m5.conf_int().loc[t, 1], 4),
                            float(f'{m5.pvalues[t]:.3g}')] for t in ['alt']}
K['R5_gencode_wide_n'] = [int(len(u)), int(u.gene.nunique()), round(u.gw_const.mean(), 3)]
K['R5_inclusion_slope'] = {t: [round(m5c.params[t], 4), float(f'{m5c.pvalues[t]:.3g}')] for t in ['incl']}

# ---------------------------------------------------------------- R6 exon-boundary placement null
ex = A.sort_values(['iso', 'exon_rank'])
bnd = {i: d.aa_end.values[:-1].astype(int) for i, d in ex.groupby('iso')}
plen = ex.groupby('iso').aa_end.max().astype(int).to_dict(); isotf = ex.groupby('iso').is_tf.first().to_dict()
segs = []
for f_, kind, s_, e_ in [('idr_exon_overlap_table.parquet', 'IDR', 'idr_aa_start', 'idr_aa_end'),
                         ('nonidr_exon_overlap_table.parquet', 'ordered', 'nonidr_aa_start', 'nonidr_aa_end')]:
    t = pd.read_parquet(os.path.join(GEN, 'interim', f_))[['isoform_accession', s_, e_]]
    t.columns = ['iso', 's', 'e']; t['kind'] = kind; segs.append(t)
S = pd.concat(segs, ignore_index=True); S = S[S.iso.isin(bnd)]; S['L'] = S.e - S.s + 1
S = S[S.L >= 20].sort_values(['iso', 's']).reset_index(drop=True)
obs, exp_ = [], []
for iso, d in S.groupby('iso', sort=True):
    Bd, P = bnd[iso], plen[iso]
    obs += list(((Bd[None, :] >= d.s.values[:, None]) & (Bd[None, :] < d.e.values[:, None])).sum(1))
    for L in d.L.values:
        n = P - L + 1
        if n <= 0: exp_.append(np.nan); continue
        st = np.arange(1, n + 1)
        exp_.append((np.searchsorted(Bd, st + L - 1, 'left') - np.searchsorted(Bd, st, 'left')).mean())
S['obs'] = obs; S['exp'] = exp_; S['tf'] = S.iso.map(isotf); S = S.dropna(subset=['exp'])
S['lbin'] = pd.cut(S.L, [19, 50, 100, 200, 400, 10 ** 5], labels=['20-50', '51-100', '101-200', '201-400', '>400'])
S.to_csv(os.path.join(OUT, 'R6_segments_obs_exp.csv'), index=False)
r6 = S.groupby(['kind', 'lbin'], observed=True).agg(n=('obs', 'size'), obs=('obs', 'sum'), exp=('exp', 'sum'),
                                                    single=('obs', lambda x: (x == 0).mean()))
r6['O_E'] = r6.obs / r6.exp; r6.to_csv(os.path.join(OUT, 'R6_by_length.csv'))
K['R6_by_length'] = {f'{a}|{b}': dict(n=int(v.n), single=round(v.single, 3), O_E=round(v.O_E, 3)) for (a, b), v in r6.iterrows()}
r6t = S.groupby('kind').agg(n=('obs', 'size'), obs=('obs', 'sum'), exp=('exp', 'sum'))
K['R6_overall_OE'] = {k: round(v.obs / v.exp, 3) for k, v in r6t.iterrows()}
pi = S.groupby(['iso', 'kind']).agg(o=('obs', 'sum'), e=('exp', 'sum')).unstack().dropna()
pi = pi[(pi[('e', 'IDR')] > 0) & (pi[('e', 'ordered')] > 0)]
ri, ro = pi[('o', 'IDR')] / pi[('e', 'IDR')], pi[('o', 'ordered')] / pi[('e', 'ordered')]
K['R6_paired_isoforms'] = dict(n=int(len(pi)), median_OE_IDR=round(ri.median(), 3), median_OE_ordered=round(ro.median(), 3),
                               wilcoxon_p=float(f'{stats.wilcoxon(ri, ro).pvalue:.3g}'))

with open(os.path.join(OUT, 'key_numbers.json'), 'w') as fh:
    json.dump(K, fh, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps(K, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
