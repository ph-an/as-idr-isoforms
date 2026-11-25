"""Disorder of every isoform (canonical and alternative) through gene age, by protein class and by protein family (Fig. 6a,b;
Extended Data Fig. 8).

For every cohort isoform: %IDR (residues in IDRs of at least 20 residues / length), the length of its longest IDR and its
mean IDR length (isoforms without an IDR have longest = 0 and no mean). Genes take their GenOrigin age and protein class.
Summaries per class (and per family with at least 15 dated genes) and age bin: number of isoforms and genes, median and
interquartile range over isoforms. Trend: Spearman correlation of gene age with the gene's mean over its isoforms (one value
per gene, so genes with many isoforms do not dominate).

Writes outputs/isoform_disorder.csv (one row per isoform), outputs/isoform_disorder_by_age.csv and
outputs/isoform_disorder_trend.csv.
"""
import os
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); OUT = os.path.join(HERE, 'outputs')
iso = pd.read_csv(os.path.join(ROOT, 'data/annotated/IDRisoforms_df.csv'), usecols=['isoform_accession', 'base_accession', 'is_canonical', 'Sequence'])
idr = pd.read_csv(os.path.join(ROOT, 'data/annotated/idr_df.csv'), usecols=['isoform_accession', 'idr_aaLen']); idr = idr[idr.idr_aaLen >= 20]
G = pd.read_csv(os.path.join(ROOT, 'analysis/protein_classes/outputs/gene_annotations.csv'), index_col=0, usecols=['gene', 'gene_age', 'age_bin', 'pclass', 'family'])
a = idr.groupby('isoform_accession').idr_aaLen.agg(total='sum', longest='max', mean_len='mean', n_idr='size')
X = iso.join(a, on='isoform_accession'); X['length'] = X.Sequence.str.len()
X[['total', 'longest', 'n_idr']] = X[['total', 'longest', 'n_idr']].fillna(0)
X['pct_idr'] = 100 * X.total / X.length
X = X.join(G, on='base_accession').dropna(subset=['gene_age'])
X.drop(columns='Sequence').to_csv(os.path.join(OUT, 'isoform_disorder.csv'), index=False)
AGE = ['<100', '100-500', '500-1000', '>1000']
FAM = X.groupby('family').base_accession.nunique(); FAM = FAM[FAM >= 15].index
rows, trend = [], []
groups = [('class', c, X[X.pclass == c]) for c in sorted(X.pclass.unique())] + [('all', 'All', X)] + [('family', f, X[X.family == f]) for f in FAM]
for kind, name, z in groups:
    for b in AGE:
        y = z[z.age_bin == b]
        r = dict(kind=kind, group=name, age_bin=b, isoforms=len(y), genes=y.base_accession.nunique())
        for m in ['pct_idr', 'longest', 'mean_len']:
            v = y[m].dropna()
            r.update({f'{m}_median': v.median() if len(v) else np.nan, f'{m}_q1': v.quantile(.25) if len(v) else np.nan, f'{m}_q3': v.quantile(.75) if len(v) else np.nan})
        rows.append(r)
    g = z.groupby('base_accession').agg(age=('gene_age', 'first'), pct_idr=('pct_idr', 'mean'), longest=('longest', 'mean'), mean_len=('mean_len', 'mean'))
    t = dict(kind=kind, group=name, genes=len(g), isoforms=len(z), age_bins=int((z.groupby('age_bin').base_accession.nunique() >= 5).sum()))
    for m in ['pct_idr', 'longest', 'mean_len']:
        gg = g.dropna(subset=[m])
        rho, p = stats.spearmanr(gg.age, gg[m]) if len(gg) > 5 else (np.nan, np.nan)
        t[f'{m}_rho'], t[f'{m}_p'] = rho, p
    trend.append(t)
pd.DataFrame(rows).to_csv(os.path.join(OUT, 'isoform_disorder_by_age.csv'), index=False)
T = pd.DataFrame(trend); T.to_csv(os.path.join(OUT, 'isoform_disorder_trend.csv'), index=False)
print(T[T.kind != 'family'].round(3).to_string()); print(T[T.kind == 'family'].sort_values('genes', ascending=False).head(20).round(3).to_string())
