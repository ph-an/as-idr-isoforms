"""With the argument `all`, every cohort isoform not already scored is predicted into outputs/v1_all/.

MetaPredict V1 (legacy network trained on MobiDB consensus disorder from sequence-based predictors, not on AlphaFold2
pLDDT) for every isoform of genes that have at least one internal exon with a lineage-restricted splice site, plus an
equal number of randomly chosen other genes. Used by sensitivity.py to ask whether the exon-age and exon-usage results
depend on a predictor trained on pLDDT. Writes outputs/v1_scores.pkl (restartable by chunk in outputs/v1/)."""
import os, pickle
import numpy as np, pandas as pd, metapredict as meta
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); OUT = os.path.join(HERE, 'outputs', 'v1'); os.makedirs(OUT, exist_ok=True)
C = pd.read_csv(os.path.join(ROOT, 'analysis/exon_age/outputs/exon_conservation.csv'), usecols=['gene', 'cons', 'pos', 'n_sites'])
young = set(C[(C.pos == 'internal') & (C.n_sites == 2) & (C.cons < 100)].gene)
others = sorted(set(C.gene) - young); rng = np.random.default_rng(20261105)
genes = young | set(rng.choice(others, len(young), replace=False))
iso = pd.read_csv(os.path.join(ROOT, 'data/annotated/IDRisoforms_df.csv'), usecols=['isoform_accession', 'base_accession', 'Sequence'])
import sys, glob
if sys.argv[1:] == ['all']:
    done = set()
    for f in glob.glob(os.path.join(OUT, 'part_*.pkl')): done |= set(pickle.load(open(f, 'rb')))
    OUT = os.path.join(HERE, 'outputs', 'v1_all'); os.makedirs(OUT, exist_ok=True)
    iso = iso[~iso.isoform_accession.isin(done)].sort_values('isoform_accession')
else:
    iso = iso[iso.base_accession.isin(genes)].sort_values('isoform_accession')
acc = iso.isoform_accession.tolist(); sq = [s.translate(str.maketrans('UBZX', 'CNQA')) for s in iso.Sequence]
print(len(young), 'young-exon genes,', len(genes), 'genes,', len(acc), 'isoforms', flush=True)
for k in range(0, len(acc), 250):
    f = os.path.join(OUT, f'part_{k // 250:03d}.pkl')
    if os.path.exists(f): continue
    r = meta.predict_disorder(dict(zip(acc[k:k + 250], sq[k:k + 250])), version='V1')
    pickle.dump({a: np.asarray(v[1] if isinstance(v, (tuple, list)) else v, dtype=np.float32) for a, v in r.items()}, open(f, 'wb')); print(k, flush=True)
