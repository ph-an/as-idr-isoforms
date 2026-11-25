"""Check that the installed MetaPredict (default settings) reproduces the project IDRs (>= 20 aa).
Re-predicts IDRs for a random sample of cohort isoforms with the installed metapredict (default disorder-domain call)
and compares segment boundaries with data/raw/metapredict/IDRome_AllSwissProtHumansProIsos.csv."""
import os, sys, numpy as np, pandas as pd, metapredict as meta
ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
iso = pd.read_csv(os.path.join(ROOT, 'data', 'annotated', 'IDRisoforms_df.csv'), usecols=['isoform_accession', 'Sequence'])
# idr_df.csv holds the project's IDRs (MetaPredict segments >= 20 aa, parsed from the precomputed IDRome per isoform)
raw = pd.read_csv(os.path.join(ROOT, 'data', 'annotated', 'idr_df.csv'), usecols=['isoform_accession', 'idr_start_aa', 'idr_end_aa'])
seg = raw.groupby('isoform_accession').apply(lambda d: sorted(zip(d.idr_start_aa, d.idr_end_aa)), include_groups=False).to_dict()
iso = iso[~iso.Sequence.str.contains('[^ACDEFGHIKLMNPQRSTVWY]', regex=True)]
s = iso.sample(N, random_state=1)
exact = res_agree = res_total = 0
for a, q in zip(s.isoform_accession, s.Sequence):
    pred = [(b[0] + 1, b[1]) for b in meta.predict_disorder_domains(q).disordered_domain_boundaries if b[1] - b[0] >= 20]
    ref = [tuple(int(x) for x in t) for t in seg.get(a, [])]
    exact += pred == ref
    m1 = np.zeros(len(q) + 1, bool); m2 = np.zeros(len(q) + 1, bool)
    for x, y in pred: m1[x:y + 1] = True
    for x, y in ref: m2[x:y + 1] = True
    res_agree += (m1[1:] == m2[1:]).sum(); res_total += len(q)
print(f'metapredict {meta.__version__}: identical segment sets for {exact}/{N} isoforms; per-residue agreement {res_agree/res_total:.4f}')
