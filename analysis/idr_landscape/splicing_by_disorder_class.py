"""Splicing-disorder coupling across the disorder spectrum: very low (<=15% IDR), intermediate, high (>=75%) genes.

Reads the outputs of idr_landscape.py (gene classes), ../reanalysis_2026-10-26 (unique exons) and
../idr_architecture (canonical-alternative pairs). Writes outputs/splicing_by_class.json.
"""
import os, json
import numpy as np, pandas as pd
from scipy import stats
H = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(H, '..')
cls = pd.read_csv(os.path.join(H, 'outputs', 'gene_disorder_class.csv'), index_col=0).dclass
ux = pd.read_csv(os.path.join(A, 'reanalysis_2026-10-26', 'outputs', 'unique_exons.csv'))
P = pd.read_csv(os.path.join(A, 'idr_architecture', 'outputs', 'pairs.csv'))
F = pd.read_csv(os.path.join(A, 'idr_architecture', 'outputs', 'canonical_idr_fate.csv'))
rng = np.random.default_rng(20261026); K = {}
ux['dclass'] = ux.gene.map(cls)
idx = ux.groupby('gene').indices
def diff(rows):
    z = ux.iloc[rows]; return z[z.alt == 1].idr_frac.mean() - z[z.alt == 0].idr_frac.mean()
K['exon_alt_minus_const'] = {}
for c in ['very low', 'intermediate', 'high']:
    genes = ux[ux.dclass == c].gene.unique()
    pt = diff(np.concatenate([idx[g] for g in genes]))
    bs = [diff(np.concatenate([idx[g] for g in rng.choice(genes, len(genes))])) for _ in range(1000)]
    z = ux[ux.dclass == c]
    K['exon_alt_minus_const'][c] = dict(genes=int(len(genes)), diff=round(pt, 4), lo=round(np.percentile(bs, 2.5), 4), hi=round(np.percentile(bs, 97.5), 4),
                                         mean_alt=round(z[z.alt == 1].idr_frac.mean(), 4), mean_const=round(z[z.alt == 0].idr_frac.mean(), 4))
# relative scale: enrichment of disorder in alternative exons relative to the gene's own constitutive exons (log2 ratio)
P['dclass'] = P.gene.map(cls)
K['spliced_sequence_disorder'] = {}
for c in ['very low', 'intermediate', 'high']:
    s = P[P.dclass == c]
    tot = lambda col, n: float((s[col] * s[n]).sum() / s[n].sum())
    K['spliced_sequence_disorder'][c] = dict(pairs=len(s), spliced_out=round(tot('dis_c_only', 'c_only'), 4), spliced_in=round(tot('dis_a_only', 'a_only'), 4),
                                             shared=round(tot('dis_shared_c', 'shared'), 4))
F['dclass'] = F.gene.map(cls)
K['canonical_idr_fate_by_class_pct'] = (pd.crosstab(F.dclass, F.fate, normalize='index') * 100).round(1).to_dict(orient='index')
json.dump(K, open(os.path.join(H, 'outputs', 'splicing_by_class.json'), 'w'), indent=1)
print(json.dumps(K, indent=1))
