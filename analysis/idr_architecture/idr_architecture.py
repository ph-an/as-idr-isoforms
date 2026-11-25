"""How alternative splicing reorganizes IDR architecture and sequence: protein-level comparison of every alternative
isoform with its gene's canonical isoform.

No transcript mapping is needed, so every isoform of the cohort enters. For each canonical-alternative pair the two
sequences are aligned (difflib, no junk heuristic) and each residue is classified as shared, canonical-only (spliced out
relative to the canonical) or alternative-only (spliced in). IDRs are MetaPredict segments >= 20 aa (data/annotated).

Outputs (outputs/):
  pairs.csv              one row per pair: residue accounting, IDR length/count change, disorder of spliced vs shared
  canonical_idr_fate.csv one row per (pair, canonical IDR): fraction retained and fate class
  alt_idr_origin.csv     one row per (pair, alternative IDR): fraction mapped to the canonical and origin class
  idr_properties_by_fate.csv  sequence properties of whole canonical IDRs by fate, with gene-clustered models
  key_numbers.json
Fate classes for a canonical IDR (r = fraction of its residues present in the alternative isoform):
  intact r = 1; removed r = 0; truncated 0 < r < 1 with the lost residues at one end of the IDR; internally deleted
  0 < r < 1 otherwise. An alternative IDR with no residue shared with the canonical is "new"; with some, "remodelled"
  if its length differs from the overlapping canonical disorder by >= 5 residues, else "carried over".
"""
import os, json, difflib
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats

ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
ANN = os.path.join(ROOT, 'data', 'annotated')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs'); os.makedirs(OUT, exist_ok=True)
iso = pd.read_csv(os.path.join(ANN, 'IDRisoforms_df.csv'), usecols=['isoform_accession', 'base_accession', 'is_canonical', 'Sequence', 'pct_idr', 'is_tf'])
idr = pd.read_csv(os.path.join(ANN, 'idr_df.csv'))
idr = idr[idr.idr_aaLen >= 20]
segs = {a: list(zip(d.idr_start_aa, d.idr_end_aa)) for a, d in idr.groupby('isoform_accession')}
seq = dict(zip(iso.isoform_accession, iso.Sequence))

def mask(acc):
    m = np.zeros(len(seq[acc]), bool)
    for a, b in segs.get(acc, []): m[a - 1:b] = True
    return m

def runs(x):
    """start/end (inclusive, 0-based) of True runs."""
    d = np.diff(np.concatenate([[0], x.astype(int), [0]]))
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0] - 1))

AA = {'aro': 'FWY', 'acid': 'DE', 'basic': 'KR', 'pro': 'P', 'st': 'ST', 'qn': 'QN', 'gly': 'G', 'hyd': 'AILMV', 'his': 'H', 'cys': 'C'}
def comp(s):
    n = len(s); r = {k: sum(s.count(c) for c in v) / n for k, v in AA.items()}
    r['fcr'] = r['acid'] + r['basic']; r['ncpr'] = r['basic'] - r['acid']; return r

pairs, fates, origins, chunks = [], [], [], []
for gene, d in iso.groupby('base_accession'):
    can = d[d.is_canonical]
    if can.empty: continue
    c = can.isoform_accession.iloc[0]; cs = seq[c]; cm = mask(c)
    for a in d[~d.is_canonical].isoform_accession:
        as_ = seq[a]; am = mask(a)
        keepC = np.zeros(len(cs), bool); keepA = np.zeros(len(as_), bool); mapC = np.full(len(cs), -1)
        for blk in difflib.SequenceMatcher(None, cs, as_, autojunk=False).get_matching_blocks():
            keepC[blk.a:blk.a + blk.size] = True; keepA[blk.b:blk.b + blk.size] = True
            mapC[blk.a:blk.a + blk.size] = np.arange(blk.b, blk.b + blk.size)
        shared = keepC.sum()
        flipsC = (cm[keepC] != am[mapC[keepC]]).sum() if shared else 0
        pairs.append(dict(gene=gene, canonical=c, alternative=a, len_c=len(cs), len_a=len(as_), shared=int(shared),
                          c_only=int((~keepC).sum()), a_only=int((~keepA).sum()),
                          dis_c_only=float(cm[~keepC].mean()) if (~keepC).any() else np.nan,
                          dis_shared_c=float(cm[keepC].mean()) if shared else np.nan,
                          dis_a_only=float(am[~keepA].mean()) if (~keepA).any() else np.nan,
                          dis_shared_a=float(am[keepA].mean()) if keepA.any() else np.nan,
                          idr_res_c=int(cm.sum()), idr_res_a=int(am.sum()), n_idr_c=len(segs.get(c, [])), n_idr_a=len(segs.get(a, [])),
                          shared_status_flips=int(flipsC), is_tf=bool(d.is_tf.iloc[0])))
        for (s0, e0) in segs.get(c, []):
            r = keepC[s0 - 1:e0]; f = r.mean()
            if f == 1: cls = 'intact'
            elif f == 0: cls = 'removed'
            else:
                lost = runs(~r); L = len(r)
                cls = 'truncated' if all(st == 0 or en == L - 1 for st, en in lost) and len(lost) <= 2 else 'internal deletion'
            pos = ((s0 + e0) / 2) / len(cs)
            fates.append(dict(gene=gene, canonical=c, alternative=a, idr_start=s0, idr_end=e0, idr_len=e0 - s0 + 1,
                              frac_retained=f, fate=cls, rel_pos=pos))
        for (s0, e0) in segs.get(a, []):
            r = keepA[s0 - 1:e0]; f = r.mean()
            if f == 0: cls = 'new'
            else:
                # canonical disorder overlapping the mapped residues of this alternative IDR
                inv = np.full(len(as_), -1); inv[mapC[keepC]] = np.where(keepC)[0]
                cpos = inv[s0 - 1:e0]; cpos = cpos[cpos >= 0]
                cov = 0
                for cs0, ce0 in runs(cm):
                    if ((cpos >= cs0) & (cpos <= ce0)).any(): cov += ce0 - cs0 + 1
                cls = 'carried over' if abs((e0 - s0 + 1) - cov) < 5 else 'remodelled'
            origins.append(dict(gene=gene, canonical=c, alternative=a, idr_len=e0 - s0 + 1, frac_shared=f, origin=cls))
        # disordered spliced vs shared chunks for composition (>= 10 residues of each kind per pair)
        for label, s_, m_, k_ in [('canonical-only', cs, cm, ~keepC), ('alternative-only', as_, am, ~keepA)]:
            sel = m_ & k_
            if sel.sum() >= 10: chunks.append(dict(gene=gene, pair=a, kind=label, where='all', n=int(sel.sum()), **comp(''.join(np.array(list(s_))[sel]))))
            # split by where the unshared run sits: at a protein terminus (new/lost N- or C-terminus) or internal
            term = np.zeros(len(s_), bool)
            for st, en in runs(k_):
                if st == 0 or en == len(s_) - 1: term[st:en + 1] = True
            for where, w in [('terminal', term), ('internal', ~term)]:
                sel2 = sel & w
                if sel2.sum() >= 10: chunks.append(dict(gene=gene, pair=a, kind=label, where=where, n=int(sel2.sum()), **comp(''.join(np.array(list(s_))[sel2]))))
        sh = cm & keepC
        if sh.sum() >= 10: chunks.append(dict(gene=gene, pair=a, kind='shared', where='all', n=int(sh.sum()), **comp(''.join(np.array(list(cs))[sh]))))
        # residue accounting by location for the alternative-only sequence
        termA = np.zeros(len(as_), bool)
        for st, en in runs(~keepA):
            if st == 0 or en == len(as_) - 1: termA[st:en + 1] = True
        pairs[-1].update(a_only_terminal=int((~keepA & termA).sum()), a_only_internal=int((~keepA & ~termA).sum()),
                         dis_a_only_terminal=float(am[~keepA & termA].mean()) if (~keepA & termA).any() else np.nan,
                         dis_a_only_internal=float(am[~keepA & ~termA].mean()) if (~keepA & ~termA).any() else np.nan)
        termC = np.zeros(len(cs), bool)
        for st, en in runs(~keepC):
            if st == 0 or en == len(cs) - 1: termC[st:en + 1] = True
        pairs[-1].update(c_only_terminal=int((~keepC & termC).sum()), c_only_internal=int((~keepC & ~termC).sum()),
                         dis_c_only_terminal=float(cm[~keepC & termC].mean()) if (~keepC & termC).any() else np.nan,
                         dis_c_only_internal=float(cm[~keepC & ~termC].mean()) if (~keepC & ~termC).any() else np.nan)

P = pd.DataFrame(pairs); F = pd.DataFrame(fates); O = pd.DataFrame(origins); Ch = pd.DataFrame(chunks)
P.to_csv(os.path.join(OUT, 'pairs.csv'), index=False); F.to_csv(os.path.join(OUT, 'canonical_idr_fate.csv'), index=False)
O.to_csv(os.path.join(OUT, 'alt_idr_origin.csv'), index=False); Ch.to_csv(os.path.join(OUT, 'idr_chunk_composition.csv'), index=False)
K = dict(n_genes=int(P.gene.nunique()), n_pairs=len(P), n_canonical_idr_involvements=len(F), n_alt_idrs=len(O))

# 1. Is spliced sequence more disordered than shared sequence? (per gene, then paired)
g = P.groupby('gene').agg(dis_c_only=('dis_c_only', 'mean'), dis_shared_c=('dis_shared_c', 'mean'),
                          dis_a_only=('dis_a_only', 'mean'), dis_shared_a=('dis_shared_a', 'mean'))
for a_, b_, lab in [('dis_c_only', 'dis_shared_c', 'spliced_out_vs_shared'), ('dis_a_only', 'dis_shared_a', 'spliced_in_vs_shared')]:
    x = g[[a_, b_]].dropna()
    K[lab] = dict(n_genes=len(x), mean_spliced=round(x[a_].mean(), 4), mean_shared=round(x[b_].mean(), 4),
                  frac_genes_higher=round((x[a_] > x[b_]).mean(), 3), wilcoxon_p=float(f'{stats.wilcoxon(x[a_], x[b_]).pvalue:.3g}'))
# residue-pooled for reference
pw = lambda col, n: round(float((P[col] * P[n]).sum() / P[n].sum()), 4)
K['pooled_fraction_disordered_by_location'] = dict(a_only_terminal=pw('dis_a_only_terminal', 'a_only_terminal'), a_only_internal=pw('dis_a_only_internal', 'a_only_internal'),
    c_only_terminal=pw('dis_c_only_terminal', 'c_only_terminal'), c_only_internal=pw('dis_c_only_internal', 'c_only_internal'),
    a_only_terminal_residue_share=round(P.a_only_terminal.sum() / P.a_only.sum(), 3), c_only_terminal_residue_share=round(P.c_only_terminal.sum() / P.c_only.sum(), 3))
K['pooled_fraction_disordered'] = dict(c_only=round((P.dis_c_only * P.c_only).sum() / P.c_only.sum(), 4),
                                       a_only=round((P.dis_a_only * P.a_only).sum() / P.a_only.sum(), 4),
                                       shared=round((P.dis_shared_c * P.shared).sum() / P.shared.sum(), 4))
# 2. IDR content change per pair
P['d_idr_res'] = P.idr_res_a - P.idr_res_c; P['d_n_idr'] = P.n_idr_a - P.n_idr_c; P['d_len'] = P.len_a - P.len_c
K['pair_changes'] = dict(frac_pairs_idr_change=round((P.d_idr_res != 0).mean(), 3), median_d_idr_res=float(P.d_idr_res.median()),
                         frac_lose_idr_res=round((P.d_idr_res < 0).mean(), 3), frac_gain_idr_res=round((P.d_idr_res > 0).mean(), 3),
                         frac_change_n_idr=round((P.d_n_idr != 0).mean(), 3),
                         frac_pairs_with_status_flip=round((P.shared_status_flips > 0).mean(), 3),
                         flipped_shared_residues_pct=round(100 * P.shared_status_flips.sum() / P.shared.sum(), 3))
# 3. Fate of canonical IDRs and origin of alternative IDRs
K['canonical_idr_fate_pct'] = (F.fate.value_counts(normalize=True) * 100).round(2).to_dict()
K['alt_idr_origin_pct'] = (O.origin.value_counts(normalize=True) * 100).round(2).to_dict()
F['pos3'] = pd.cut(F.rel_pos, [0, 1 / 3, 2 / 3, 1.0001], labels=['N-terminal', 'middle', 'C-terminal'])
K['fate_by_position_pct'] = (pd.crosstab(F.pos3, F.fate, normalize='index') * 100).round(1).to_dict(orient='index')
F['lbin'] = pd.cut(F.idr_len, [19, 50, 100, 200, 10 ** 5], labels=['20-50', '51-100', '101-200', '>200'])
K['fate_by_length_pct'] = (pd.crosstab(F.lbin, F.fate, normalize='index') * 100).round(1).to_dict(orient='index')
# 4. Sequence properties of whole canonical IDRs: removed vs intact (same canonical isoform; length-adjusted)
props = ['FCR', 'NCPR', 'kappa', 'fract_aro', 'fract_pro', 'fract_polar', 'fract_ali', 'fract_negative', 'fract_positive', 'scaling_exponent']
cid = idr[idr.isoform_accession.isin(iso[iso.is_canonical].isoform_accession)][['isoform_accession', 'idr_start_aa', 'idr_end_aa', 'idr_aaLen'] + props]
ever = F.groupby(['canonical', 'idr_start', 'idr_end']).fate.agg(lambda s: 'removed' if (s == 'removed').any() else ('intact' if (s == 'intact').all() else 'altered')).reset_index()
cid = cid.merge(ever, left_on=['isoform_accession', 'idr_start_aa', 'idr_end_aa'], right_on=['canonical', 'idr_start', 'idr_end'])
cid['loglen'] = np.log10(cid.idr_aaLen); cid['removed'] = (cid.fate == 'removed').astype(int)
sub = cid[cid.fate.isin(['removed', 'intact'])]
rows = []
for p in props:
    m = smf.ols(f'{p} ~ removed + loglen', sub).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(sub.isoform_accession)[0]})
    rows.append(dict(property=p, mean_intact=sub[sub.removed == 0][p].mean(), mean_removed=sub[sub.removed == 1][p].mean(),
                     adj_diff=m.params.removed, lo=m.conf_int().loc['removed', 0], hi=m.conf_int().loc['removed', 1], p=m.pvalues.removed))
R = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
R['q_BH'] = multipletests(R.p, method='fdr_bh')[1]
R.to_csv(os.path.join(OUT, 'idr_properties_by_fate.csv'), index=False)
K['n_canonical_idrs_removed_vs_intact'] = dict(removed=int(sub.removed.sum()), intact=int((sub.removed == 0).sum()))
K['idr_properties_removed_vs_intact'] = R.round(4).to_dict(orient='records')
# 5. Composition of spliced (alternative- or canonical-only) disorder vs shared disorder within pairs
feats = [k for k in AA] + ['fcr', 'ncpr']
crow = []
for kind, where in [(k, w) for k in ['canonical-only', 'alternative-only'] for w in ['all', 'internal', 'terminal']]:
    x = Ch[((Ch.kind == kind) & (Ch['where'] == where)) | (Ch.kind == 'shared')].copy(); x['spliced'] = (x.kind == kind).astype(int); x['logn'] = np.log10(x.n)
    for f in feats:
        m = smf.ols(f'{f} ~ spliced + logn', x).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(x.gene)[0]})
        crow.append(dict(kind=kind, where=where, feature=f, mean_shared=x[x.spliced == 0][f].mean(), mean_spliced=x[x.spliced == 1][f].mean(),
                         adj_diff_pp=100 * m.params.spliced, lo=100 * m.conf_int().loc['spliced', 0], hi=100 * m.conf_int().loc['spliced', 1], p=m.pvalues.spliced))
C = pd.DataFrame(crow); C['q_BH'] = C.groupby(['kind', 'where']).p.transform(lambda p: multipletests(p, method='fdr_bh')[1])
C.to_csv(os.path.join(OUT, 'spliced_disorder_composition.csv'), index=False)
K['spliced_disorder_composition'] = C.round(4).to_dict(orient='records')
json.dump(K, open(os.path.join(OUT, 'key_numbers.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps({k: v for k, v in K.items() if k not in ('idr_properties_removed_vs_intact', 'spliced_disorder_composition')}, indent=1, default=str))
print(R.round(4).to_string()); print(C[C.q_BH < 0.05].round(4).to_string())
