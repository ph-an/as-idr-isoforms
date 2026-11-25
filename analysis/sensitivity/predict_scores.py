"""Per-residue MetaPredict (v3.1) disorder scores for every cohort isoform, so that IDR calls can be re-derived at other
thresholds and minimum lengths (threshold_sensitivity.py). Writes outputs/scores/part_###.pkl in chunks of 500 isoforms
(restartable; existing parts are skipped)."""
import os, pickle, sys
import numpy as np, pandas as pd, metapredict as meta
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, 'outputs', 'scores'); os.makedirs(OUT, exist_ok=True)
iso = pd.read_csv(os.path.join(ROOT, 'data/annotated/IDRisoforms_df.csv'), usecols=['isoform_accession', 'Sequence']).sort_values('isoform_accession')
acc = iso.isoform_accession.tolist()
# MetaPredict accepts the 20 standard residues only; the few non-standard ones are replaced (U->C, B->N, Z->Q, X->A)
sq = [s.translate(str.maketrans('UBZX', 'CNQA')) for s in iso.Sequence]
for k in range(0, len(acc), 500):
    f = os.path.join(OUT, f'part_{k // 500:03d}.pkl')
    if os.path.exists(f): continue
    chunk = dict(zip(acc[k:k + 500], sq[k:k + 500]))
    short = {a: q for a, q in chunk.items() if len(q) <= 3000}   # long proteins one at a time (padding a batch to titin length runs out of memory)
    r = meta.predict_disorder(short, batch_size=32) if short else {}
    for a, q in chunk.items():
        if a not in short: r[a] = meta.predict_disorder(q)
    pickle.dump({a: np.asarray(v[1] if isinstance(v, (tuple, list)) else v, dtype=np.float32) for a, v in r.items()}, open(f, 'wb'))
    print(k, flush=True)
