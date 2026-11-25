"""Where does disordered sequence found only in alternative isoforms come from in the genome?

For canonical-alternative pairs in which both isoforms have a validated exon structure, each residue of the
alternative isoform is mapped to the genomic position of the first base of its codon. Residues absent from the
canonical isoform (protein alignment, as in idr_architecture.py) and internal (flanked by shared sequence) are
classified by that position:
  novel       - not inside the canonical isoform's coding sequence (alternative exon or extension)
  same frame  - inside the canonical CDS and in the same codon position (should be rare: same codon, aligned away)
  other frame - inside the canonical CDS but in a different reading frame (frame-shifted reuse of canonical exons)
The composition of disordered residues in each class is compared with the disordered residues shared by the pair
(adjusted for log pool size, gene-clustered errors). Writes outputs/spliced_sequence_origin.json.
"""
import os, json, difflib
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')
AA = {'aro': 'FWY', 'acid': 'DE', 'basic': 'KR', 'pro': 'P', 'st': 'ST', 'qn': 'QN', 'gly': 'G', 'hyd': 'AILMV', 'his': 'H', 'cys': 'C'}
E = pd.read_parquet(os.path.join(ROOT, 'data/genomic/processed/exon_level_table.parquet'),
                    columns=['isoform_accession', 'base_accession', 'is_canonical', 'is_coding_exon', 'strand', 'cds_start', 'cds_end', 'aa_start', 'chr'])
E = E[E.is_coding_exon & E.cds_start.notna()]
iso = pd.read_csv(os.path.join(ROOT, 'data/annotated/IDRisoforms_df.csv'), usecols=['isoform_accession', 'Sequence'])
seq = dict(zip(iso.isoform_accession, iso.Sequence))
idr = pd.read_csv(os.path.join(ROOT, 'data/annotated/idr_df.csv')); idr = idr[idr.idr_aaLen >= 20]
segs = {a: list(zip(d.idr_start_aa, d.idr_end_aa)) for a, d in idr.groupby('isoform_accession')}
def mask(a):
    m = np.zeros(len(seq[a]), bool)
    for s, e in segs.get(a, []): m[int(s) - 1:int(e)] = True
    return m
def coding_bases(d):
    d = d.sort_values('aa_start', kind='stable'); strand = d.strand.iloc[0]
    parts = [np.arange(int(r.cds_start), int(r.cds_end) + 1) if strand == '+' else np.arange(int(r.cds_end), int(r.cds_start) - 1, -1) for r in d.itertuples()]
    return np.concatenate(parts)
def comp(s):
    n = len(s); c = {k: sum(s.count(x) for x in v) / n for k, v in AA.items()}
    c['fcr'] = c['acid'] + c['basic']; c['ncpr'] = (c['basic'] - c['acid']); return c

rows, chunks = [], []
for gene, d in E.groupby('base_accession'):
    if not d.is_canonical.any(): continue
    c = d[d.is_canonical].isoform_accession.iloc[0]
    if c not in seq: continue
    cb = coding_bases(d[d.isoform_accession == c])
    if len(cb) != 3 * len(seq[c]): continue
    cpos = {g: i % 3 for i, g in enumerate(cb)}  # genomic base -> codon position in canonical
    cm = mask(c)
    for a in d[~d.is_canonical].isoform_accession.unique():
        if a not in seq: continue
        ab = coding_bases(d[d.isoform_accession == a])
        if len(ab) != 3 * len(seq[a]): continue
        am = mask(a); keepA = np.zeros(len(seq[a]), bool); keepC = np.zeros(len(seq[c]), bool)
        for blk in difflib.SequenceMatcher(None, seq[c], seq[a], autojunk=False).get_matching_blocks():
            keepA[blk.b:blk.b + blk.size] = True; keepC[blk.a:blk.a + blk.size] = True
        if not keepA.any(): continue
        first, last = np.argmax(keepA), len(keepA) - 1 - np.argmax(keepA[::-1])
        internal = np.zeros(len(keepA), bool); internal[first:last + 1] = True
        sel = ~keepA & internal
        if not sel.any(): continue
        cls = np.array(['novel' if (g := ab[3 * i]) not in cpos else ('same frame' if cpos[g] == 0 else 'other frame') for i in np.where(sel)[0]])
        idxs = np.where(sel)[0]
        for k in ['novel', 'same frame', 'other frame']:
            r = idxs[cls == k]
            rows.append(dict(gene=gene, pair=a, cls=k, n=len(r), n_dis=int(am[r].sum())))
            rd = r[am[r]]
            if len(rd) >= 10: chunks.append(dict(gene=gene, pair=a, kind=k, n=len(rd), **comp(''.join(seq[a][i] for i in rd))))
        sh = keepC & cm
        if sh.sum() >= 10: chunks.append(dict(gene=gene, pair=a, kind='shared', n=int(sh.sum()), **comp(''.join(np.array(list(seq[c]))[sh]))))
R = pd.DataFrame(rows); Ch = pd.DataFrame(chunks)
K = {'pairs': int(R.pair.nunique()), 'genes': int(R.gene.nunique())}
t = R.groupby('cls')[['n', 'n_dis']].sum()
K['internal_alt_only_residues'] = {k: dict(residues=int(v.n), disordered=int(v.n_dis), frac_disordered=round(v.n_dis / v.n, 4) if v.n else None,
                                           share_of_disordered=round(v.n_dis / t.n_dis.sum(), 4)) for k, v in t.iterrows()}
feats = list(AA) + ['fcr', 'ncpr']; crow = []
for k in ['novel', 'other frame']:
    x = Ch[Ch.kind.isin([k, 'shared'])].copy(); x['s'] = (x.kind == k).astype(int); x['logn'] = np.log10(x.n)
    if x.s.sum() < 30: continue
    for f in feats:
        m = smf.ols(f'{f} ~ s + logn', x).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(x.gene)[0]})
        crow.append(dict(kind=k, feature=f, n=int(x.s.sum()), diff_pp=100 * m.params.s, lo=100 * m.conf_int().loc['s', 0], hi=100 * m.conf_int().loc['s', 1], p=m.pvalues.s))
C = pd.DataFrame(crow); C['q_BH'] = C.groupby('kind').p.transform(lambda p: multipletests(p, method='fdr_bh')[1])
C.to_csv(os.path.join(OUT, 'spliced_sequence_origin_composition.csv'), index=False)
K['composition'] = C.round(4).to_dict(orient='records')
json.dump(K, open(os.path.join(OUT, 'spliced_sequence_origin.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in K.items() if k != 'composition'}, indent=1)); print(C.round(3).to_string())
