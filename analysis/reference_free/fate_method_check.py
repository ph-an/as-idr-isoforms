"""Do IDR fates depend on how 'shared' residues are defined, or on the reference isoform?

Compares, for the same canonical IDR in the same canonical-alternative pair, the fate assigned by protein-sequence
alignment (analysis/idr_architecture, difflib) with the fate assigned by shared genomic codons (the definition used in
reference_free.py), and computes genomic-codon fates in both directions (canonical->alternative, alternative->canonical).
Writes outputs/fate_method_check.json and outputs/canonical_fate_alignment_vs_genomic.csv.
"""
import os, json, collections
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(HERE, '..', '..')))
exec(open(os.path.join(HERE, 'reference_free.py')).read().split('gene_rows, chunk_rows')[0])

def classify(x):
    f = x.mean()
    if f == 1: return 'intact'
    if f == 0: return 'removed'
    lost = runs(~x); L = len(x)
    return 'truncated' if all(st == 0 or en == L - 1 for st, en in lost) and len(lost) <= 2 else 'internal deletion'

keys = {}
for a, da in E.groupby('isoform_accession'):
    if a in seq:
        k, nb = codon_keys(da)
        if nb == 3 * len(seq[a]): keys[a] = k
sets = {a: set(k.tolist()) for a, k in keys.items()}
K = {}
Fc = pd.read_csv(os.path.join(HERE, '..', 'idr_architecture', 'outputs', 'canonical_idr_fate.csv'))
Fc = Fc[Fc.canonical.isin(keys) & Fc.alternative.isin(keys)].copy()
Fc['genomic_fate'] = [classify(np.fromiter((v in sets[r.alternative] for v in keys[r.canonical][r.idr_start - 1:r.idr_end]), bool)) for r in Fc.itertuples()]
Fc.to_csv(os.path.join(OUT, 'canonical_fate_alignment_vs_genomic.csv'), index=False)
K['same_idrs'] = len(Fc)
K['alignment_fate_pct'] = (100 * Fc.fate.value_counts(normalize=True)).round(1).to_dict()
K['genomic_fate_pct'] = (100 * Fc.genomic_fate.value_counts(normalize=True)).round(1).to_dict()
K['crosstab_alignment_rows_genomic_cols'] = pd.crosstab(Fc.fate, Fc.genomic_fate).to_dict(orient='index')
sub = Fc[(Fc.fate == 'internal deletion') & (Fc.genomic_fate != 'internal deletion')]
K['alignment_internal_deletions_not_genomic'] = dict(n=len(sub), of=int((Fc.fate == 'internal deletion').sum()), genomic=sub.genomic_fate.value_counts().to_dict())
cnt = collections.Counter()
for gene, d in E.groupby('base_accession'):
    can = [a for a in d[d.is_canonical].isoform_accession.unique() if a in keys]
    alts = [a for a in d[~d.is_canonical].isoform_accession.unique() if a in keys]
    if not can: continue
    c = can[0]
    for a in alts:
        for src, dst, lab in [(c, a, 'canonical_to_alternative'), (a, c, 'alternative_to_canonical')]:
            for s, e in segs.get(src, []):
                cnt[(lab, classify(np.fromiter((v in sets[dst] for v in keys[src][s - 1:e]), bool)))] += 1
D = pd.Series(cnt).unstack()
K['genomic_fate_pct_by_direction'] = (100 * D.div(D.sum(axis=1), axis=0)).round(1).to_dict(orient='index')
json.dump(K, open(os.path.join(OUT, 'fate_method_check.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps(K, indent=1, default=str))
