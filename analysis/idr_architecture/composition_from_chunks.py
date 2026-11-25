"""Recompute section 5 of idr_architecture.py (composition of spliced vs shared disorder) from the saved
outputs/idr_chunk_composition.csv, without re-running the alignments, and update key_numbers.json."""
import os, json
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')
Ch = pd.read_csv(os.path.join(OUT, 'idr_chunk_composition.csv'))
feats = ['aro', 'acid', 'basic', 'pro', 'st', 'qn', 'gly', 'hyd', 'his', 'cys', 'fcr', 'ncpr']
crow = []
for kind in ['canonical-only', 'alternative-only']:
    for where in ['all', 'internal', 'terminal']:
        x = Ch[((Ch.kind == kind) & (Ch['where'] == where)) | (Ch.kind == 'shared')].copy()
        x['spliced'] = (x.kind == kind).astype(int); x['logn'] = np.log10(x.n)
        for f in feats:
            m = smf.ols(f'{f} ~ spliced + logn', x).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(x.gene)[0]})
            crow.append(dict(kind=kind, where=where, feature=f, n_spliced=int(x.spliced.sum()), mean_shared=x[x.spliced == 0][f].mean(),
                             mean_spliced=x[x.spliced == 1][f].mean(), adj_diff_pp=100 * m.params.spliced,
                             lo=100 * m.conf_int().loc['spliced', 0], hi=100 * m.conf_int().loc['spliced', 1], p=m.pvalues.spliced))
C = pd.DataFrame(crow); C['q_BH'] = C.groupby(['kind', 'where']).p.transform(lambda p: multipletests(p, method='fdr_bh')[1])
C.to_csv(os.path.join(OUT, 'spliced_disorder_composition.csv'), index=False)
K = json.load(open(os.path.join(OUT, 'key_numbers.json'))); K['spliced_disorder_composition'] = C.round(4).to_dict(orient='records')
json.dump(K, open(os.path.join(OUT, 'key_numbers.json'), 'w'), indent=1)
print(C[C.q_BH < 0.05].round(3).to_string())

# Sensitivity: alternative-only internal disorder vs shared disorder, excluding genes with UniProt keyword 'Zinc-finger'
ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
kw = pd.read_csv(os.path.join(ROOT, 'data', 'raw', 'uniprot', 'uniprotkb_annotated'), sep='\t', usecols=['Entry', 'Keywords'], low_memory=False).set_index('Entry').Keywords
zf = set(kw.index[kw.fillna('').str.contains('Zinc-finger')])
x = Ch[(((Ch.kind == 'alternative-only') & (Ch['where'] == 'internal')) | (Ch.kind == 'shared')) & ~Ch.gene.isin(zf)].copy()
x['spliced'] = (x.kind != 'shared').astype(int); x['logn'] = np.log10(x.n)
nz = {'n_spliced': int(x.spliced.sum())}
for f in feats:
    m = smf.ols(f'{f} ~ spliced + logn', x).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(x.gene)[0]})
    nz[f] = [round(100 * m.params.spliced, 3), round(100 * m.conf_int().loc['spliced', 0], 3), round(100 * m.conf_int().loc['spliced', 1], 3), float(f'{m.pvalues.spliced:.3g}')]
K['alt_only_internal_non_zinc_finger_genes'] = nz
json.dump(K, open(os.path.join(OUT, 'key_numbers.json'), 'w'), indent=1)
print(nz)
