"""The disorder landscape of human protein isoforms (manuscript Fig. 2): isoform disorder content, IDR length vs disorder
content, protein length, and highly disordered vs barely disordered proteins.

Inputs: data/annotated/IDRisoforms_df_geneage.csv, data/annotated/idr_df.csv, data/raw/uniprot/uniprotkb_annotated.
Disorder classes use the canonical isoform's %IDR with the thresholds of Song & Akey (2026): high >= 75%,
very low <= 15%, intermediate otherwise.
"""
import os, json
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
ANN = os.path.join(ROOT, 'data', 'annotated'); RAW = os.path.join(ROOT, 'data', 'raw')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs'); os.makedirs(OUT, exist_ok=True)
HI, LO = 75, 15
d = pd.read_csv(os.path.join(ANN, 'IDRisoforms_df_geneage.csv'), low_memory=False)
d['is_tf'] = d.tf_group.eq('TF')
idr = pd.read_csv(os.path.join(ANN, 'idr_df.csv')); idr = idr[idr.idr_aaLen >= 20]
d['n_idr_segments'] = d.n_idr_segments.fillna(0); d['max_idr_len'] = d.max_idr_len.fillna(0)
K = {}
g = d.groupby('base_accession').agg(n_iso=('isoform_accession', 'size'), is_tf=('is_tf', 'first'))
K['cohort'] = dict(isoforms=len(d), genes=len(g), canonical=int(d.is_canonical.sum()), alternative=int((~d.is_canonical).sum()),
                   idrs=len(idr[idr.isoform_accession.isin(d.isoform_accession)]),
                   genes_by_isoform_count={str(k): int(v) for k, v in g.n_iso.clip(upper=6).value_counts().sort_index().items()},
                   isoforms_with_idr_pct=round(100 * (d.n_idr_segments > 0).mean(), 1),
                   residue_weighted_pct_idr=round(100 * d.idr_aaLen.fillna(0).sum() / d.Length.sum(), 2),
                   median_pct_idr=round(d.pct_idr.median(), 2))
# isoform disorder content: canonical vs alternative, paired within gene
can = d[d.is_canonical].set_index('base_accession')
alt = d[~d.is_canonical].copy(); alt['can_pct'] = alt.base_accession.map(can.pct_idr); alt['can_len'] = alt.base_accession.map(can.Length)
alt['d_pct'] = alt.pct_idr - alt.can_pct; alt['d_len'] = alt.Length - alt.can_len
K['canonical_vs_alternative'] = dict(median_canonical=round(can.pct_idr.median(), 2), median_alternative=round(alt.pct_idr.median(), 2),
    mwu_p=float(f'{stats.mannwhitneyu(can.pct_idr, alt.pct_idr).pvalue:.3g}'),
    paired_median_delta=round(alt.d_pct.median(), 2), frac_alt_lower=round((alt.d_pct < 0).mean(), 3), frac_alt_higher=round((alt.d_pct > 0).mean(), 3),
    frac_alt_equal=round((alt.d_pct == 0).mean(), 3), frac_alt_shorter=round((alt.d_len < 0).mean(), 3))
# IDR length vs disorder content; protein length vs disorder content
w = d[d.n_idr_segments > 0]
K['idr_length_vs_pct_idr'] = dict(n=len(w), spearman_longest=[round(x, 4) for x in stats.spearmanr(w.max_idr_len, w.pct_idr)],
                                  spearman_mean=[round(x, 4) for x in stats.spearmanr(w.mean_idr_len, w.pct_idr)],
                                  spearman_nsegments=[round(x, 4) for x in stats.spearmanr(w.n_idr_segments, w.pct_idr)],
                                  median_longest_idr=float(w.max_idr_len.median()), median_idr_len_all_segments=float(idr[idr.isoform_accession.isin(d.isoform_accession)].idr_aaLen.median()))
K['protein_length_vs_pct_idr'] = dict(spearman=[round(x, 4) for x in stats.spearmanr(d.Length, d.pct_idr)], median_length=float(d.Length.median()))
# disorder classes from the canonical isoform
can = can.copy(); can['dclass'] = np.select([can.pct_idr >= HI, can.pct_idr <= LO], ['high', 'very low'], 'intermediate')
gene_range = d.groupby('base_accession').pct_idr.agg(lambda s: s.max() - s.min())
can['range'] = gene_range; can['n_iso'] = g.n_iso
cls = can.groupby('dclass').agg(genes=('Length', 'size'), median_length=('Length', 'median'), median_n_iso=('n_iso', 'median'),
                                 mean_n_iso=('n_iso', 'mean'), median_range=('range', 'median'), mean_range=('range', 'mean'),
                                 frac_tf=('is_tf', 'mean'), median_age=('gene_age', 'median'))
K['disorder_classes'] = cls.round(3).to_dict(orient='index')
hl = can[can.dclass.isin(['high', 'very low'])]
K['high_vs_very_low_tests'] = {c: float(f'{stats.mannwhitneyu(hl[hl.dclass == "high"][c].dropna(), hl[hl.dclass == "very low"][c].dropna()).pvalue:.3g}')
                              for c in ['Length', 'n_iso', 'range']}
# alternative isoforms of each class: how far their disorder moves from the canonical
alt['dclass'] = alt.base_accession.map(can.dclass)
K['alt_delta_by_class'] = alt.groupby('dclass').d_pct.agg(['size', 'median', lambda s: (s.abs() >= 10).mean()]).round(3).rename(columns={'<lambda_0>': 'frac_abs_change_ge10pp'}).to_dict(orient='index')
# UniProt keyword enrichment, high and very-low classes vs the whole cohort (gene level)
kw = pd.read_csv(os.path.join(RAW, 'uniprot', 'uniprotkb_annotated'), sep='\t', usecols=['Entry', 'Keywords'], low_memory=False).set_index('Entry').Keywords
genes = can.index.intersection(kw.index); sets = {gname: set(str(kw[gname]).split(';')) if isinstance(kw[gname], str) else set() for gname in genes}
allk = pd.Series([k for s in sets.values() for k in s]).value_counts(); allk = allk[allk >= 20].index
rows = []
for c in ['high', 'very low']:
    members = set(can.index[can.dclass == c]) & set(genes); N, n = len(genes), len(members)
    for k in allk:
        K_ = sum(k in s for s in sets.values()); x = sum(k in sets[gname] for gname in members)
        o, p = stats.fisher_exact([[x, n - x], [K_ - x, N - n - (K_ - x)]], alternative='greater')
        rows.append(dict(cls=c, keyword=k, in_class=x, class_size=n, in_cohort=K_, cohort=N, fold=(x / n) / (K_ / N) if K_ else np.nan, p=p))
E = pd.DataFrame(rows); E['q_BH'] = E.groupby('cls').p.transform(lambda p: multipletests(p, method='fdr_bh')[1])
E = E.sort_values(['cls', 'q_BH']); E.to_csv(os.path.join(OUT, 'keyword_enrichment_by_disorder_class.csv'), index=False)
K['top_keywords'] = {c: E[(E.cls == c) & (E.q_BH < 0.05) & (E.in_class >= 10)].sort_values('fold', ascending=False).head(12)[['keyword', 'in_class', 'fold', 'q_BH']].round(4).to_dict(orient='records') for c in ['high', 'very low']}
can[['dclass', 'pct_idr', 'Length', 'n_iso', 'range', 'is_tf']].to_csv(os.path.join(OUT, 'gene_disorder_class.csv'))
alt[['isoform_accession', 'base_accession', 'pct_idr', 'can_pct', 'd_pct', 'd_len', 'dclass']].to_csv(os.path.join(OUT, 'alt_vs_canonical.csv'), index=False)
json.dump(K, open(os.path.join(OUT, 'key_numbers.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps(K, indent=1, default=str))
