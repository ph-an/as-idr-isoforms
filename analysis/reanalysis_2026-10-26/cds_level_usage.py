"""Sensitivity: exon usage defined on coding coordinates instead of Ensembl exon IDs.

Ensembl exon IDs encode the full exon including UTR, so two transcripts that share an identical coding segment but
differ in their UTR (alternative transcription start or polyadenylation) carry different exon IDs, and the ID-based
definition calls that coding segment "alternative". Here a coding segment is keyed by (gene, chr, CDS start, CDS end);
it is constitutive if every mapped coding transcript of the gene contains it. The main contrast (alternative vs
constitutive IDR fraction, gene bootstrap and gene-paired Wilcoxon) is recomputed on these segments.

Writes outputs/cds_level_usage.json.
"""
import os, json
import numpy as np, pandas as pd
from scipy import stats

ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')
rng = np.random.default_rng(20261026)
A = pd.read_csv(os.path.join(ROOT, 'data', 'genomic', 'analysis', 'exon_idr_annotation.csv'))
E = pd.read_parquet(os.path.join(ROOT, 'data', 'genomic', 'processed', 'exon_level_table.parquet'),
                    columns=['transcript_id', 'exon_id', 'chr', 'cds_start', 'cds_end', 'is_coding_exon'])
E = E[E.is_coding_exon & E.cds_start.notna()].drop_duplicates(['transcript_id', 'exon_id'])
A = A.merge(E.rename(columns={'transcript_id': 'ENST'}), on=['ENST', 'exon_id'], how='left')
K = {'rows': len(A), 'rows_without_cds_coords': int(A.cds_start.isna().sum())}
A = A[A.cds_start.notna() & (A.n_tx_gene >= 2)].copy()
A['seg'] = A.chr + ':' + A.cds_start.astype(int).astype(str) + '-' + A.cds_end.astype(int).astype(str)
ntx = A.groupby('gene').ENST.nunique()
seg = A.groupby(['gene', 'seg']).agg(n=('ENST', 'nunique'), idr_frac=('idr_frac', 'mean'), id_usage=('usage', lambda s: (s == 'alternative').any()),
                                     aa=('aa_end', 'size')).reset_index()
seg['ntx'] = seg.gene.map(ntx); seg = seg[seg.ntx >= 2]
seg['alt'] = (seg.n < seg.ntx).astype(int)
K['segments'] = len(seg); K['genes'] = int(seg.gene.nunique())
K['segments_alt'] = int(seg.alt.sum()); K['segments_const'] = int((1 - seg.alt).sum())
K['id_alt_but_cds_const'] = int((seg.id_usage & (seg.alt == 0)).sum())
K['mean_idr_alt'] = round(seg[seg.alt == 1].idr_frac.mean(), 4); K['mean_idr_const'] = round(seg[seg.alt == 0].idr_frac.mean(), 4)
idx = seg.groupby('gene').indices; genes = np.array(list(idx))
def diff(rows):
    z = seg.iloc[rows]; return z[z.alt == 1].idr_frac.mean() - z[z.alt == 0].idr_frac.mean()
bs = [diff(np.concatenate([idx[x] for x in rng.choice(genes, len(genes))])) for _ in range(1000)]
K['diff'] = round(K['mean_idr_alt'] - K['mean_idr_const'], 4); K['diff_ci95'] = [round(np.percentile(bs, 2.5), 4), round(np.percentile(bs, 97.5), 4)]
pg = seg.groupby(['gene', 'alt']).idr_frac.mean().unstack().dropna()
K['paired_genes'] = len(pg); K['paired_frac_alt_higher'] = round(float((pg[1] > pg[0]).mean()), 3)
K['paired_wilcoxon_p'] = float(f'{stats.wilcoxon(pg[1], pg[0]).pvalue:.3g}')
x = seg.idr_frac >= 0.5
t = pd.crosstab(seg.alt, x); K['odds_ratio_idr_exon'] = round(float(stats.contingency.odds_ratio(t.values[::-1, ::-1]).statistic), 3)
json.dump(K, open(os.path.join(OUT, 'cds_level_usage.json'), 'w'), indent=1)
print(json.dumps(K, indent=1))
