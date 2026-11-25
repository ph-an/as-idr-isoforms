"""Robustness of the IDR-fate comparison (Fig. 5b-d; Methods, reference-free comparison).

Questions
  R1  Does the alignment-vs-codon contrast depend on how protein sequences are aligned? The original alignment used
      difflib.SequenceMatcher (autojunk=False), which is a greedy longest-common-block matcher, not a biological
      alignment: every matching block, however short, counts as shared sequence. We classify the same canonical IDRs in
      the same canonical-alternative pairs with
        difflib autojunk=False (original), difflib autojunk=True, difflib keeping only blocks of >= 5 or >= 10 residues,
        global Needleman-Wunsch (Biopython PairwiseAligner, BLOSUM62, gap open -10, extend -0.5; EMBOSS needle defaults)
        counting identical aligned residues as shared, and the same alignment keeping identical runs of >= 5 residues,
      against the codon definitions below.
  R2  Junction codons. The codon key is the genomic coordinate of a codon's first base. A codon that spans an exon
      junction can share its first base with a codon of another isoform whose second and third bases come from a
      different exon. We re-key each codon by all three genomic bases, and additionally require the same amino acid,
      and count how often the first-base key and the full key disagree.
  R3  Weighting. Ordered pairs grow as n(n-1) with a gene's isoform count, so a few genes contribute many comparisons.
      Fates are summarised pooled, with every gene weighted equally, by isoform-count stratum, without the 5% of genes
      that contribute the most comparisons, and with one random canonical-alternative pair per gene; gene-bootstrap 95%
      intervals.
  R4  Does keeping an IDR's codons keep its predicted disorder? For every IDR of isoform A followed into isoform B, the
      fraction of its shared codons that are inside an IDR in B's own MetaPredict prediction.

Writes outputs/fate_robustness.json and outputs/fate_robustness_methods.csv.
"""
import os, json, difflib, itertools, time
import numpy as np, pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(HERE, '..', '..')))
exec(open(os.path.join(HERE, 'reference_free.py')).read().split('gene_rows, chunk_rows')[0])
rng = np.random.default_rng(20261105)
MAXLEN_NW = 6000  # global alignment of longer pairs is skipped (memory); the comparison table is restricted to pairs all methods cover


def classify(x):
    f = x.mean()
    if f == 1: return 'intact'
    if f == 0: return 'removed'
    lost = runs(~x); L = len(x)
    return 'truncated' if all(st == 0 or en == L - 1 for st, en in lost) and len(lost) <= 2 else 'internal deletion'


def bases(d):
    d = d.sort_values('aa_start', kind='stable'); strand = d.strand.iloc[0]
    parts = [np.arange(int(r.cds_start), int(r.cds_end) + 1) if strand == '+' else np.arange(int(r.cds_end), int(r.cds_start) - 1, -1)
             for r in d.itertuples()]
    return np.concatenate(parts)


R3_ONLY = os.environ.get('FATE_R3_ONLY') == '1'   # skip R1/R2 (reuse their saved results) and recompute R3/R4
fates = ['intact', 'truncated', 'removed', 'internal deletion']
B, gene_of, can_of = {}, {}, {}
for a, da in E.groupby('isoform_accession'):
    if a in seq:
        b = bases(da)
        if len(b) == 3 * len(seq[a]):
            B[a] = b.reshape(-1, 3).astype(np.int64); gene_of[a] = da.base_accession.iloc[0]
            if da.is_canonical.iloc[0]: can_of[gene_of[a]] = a
del E
first = {a: b[:, 0] for a, b in B.items()}
# all three genomic bases of a codon combined into one 64-bit key (a hash; keys are only compared within a gene)
with np.errstate(over='ignore'):
    full = {a: (b[:, 0].astype(np.uint64) * np.uint64(1000003) + b[:, 1].astype(np.uint64) * np.uint64(10007) + b[:, 2].astype(np.uint64)) for a, b in B.items()}
aa = {a: np.frombuffer(seq[a].encode(), dtype='S1') for a in B}
genes = pd.Series(gene_of).groupby(pd.Series(gene_of)).groups
K = {}

if R3_ONLY:
    K = json.load(open(os.path.join(OUT, 'fate_robustness.json')))
    K = {k: v for k, v in K.items() if k.startswith(('R1', 'R2'))}
else:
    # ---------------- R2: junction codons
    spl = {a: int((np.abs(np.diff(b, axis=1)) != 1).any(axis=1).sum()) for a, b in B.items()}
    K['R2_codons_spanning_a_junction'] = dict(codons=int(sum(spl.values())), of=int(sum(len(b) for b in B.values())))
    mis_key, mis_aa, shared_first = 0, 0, 0
    for gene, accs in genes.items():
        for a, b in itertools.permutations(list(accs), 2):
            common, ia, ib = np.intersect1d(first[a], first[b], return_indices=True)
            shared_first += len(common); d = full[a][ia] != full[b][ib]
            mis_key += int(d.sum()); mis_aa += int((d & (aa[a][ia] != aa[b][ib])).sum())
    K['R2_first_base_key_shared_but_other_bases_differ'] = dict(pairs_of_residues=shared_first, differ=mis_key, of_which_amino_acid_differs=mis_aa,
                                                                frac_differ=round(mis_key / shared_first, 6))
    print(K, flush=True)

    # ---------------- R1: alignment methods on the same canonical IDRs
    mat = substitution_matrices.load('BLOSUM62')
    al = Align.PairwiseAligner(); al.substitution_matrix = mat; al.open_gap_score = -10; al.extend_gap_score = -0.5; al.mode = 'global'


    def difflib_shared(s, t, autojunk, minblk=1):
        kc = np.zeros(len(s), bool); m = np.full(len(s), -1)
        for blk in difflib.SequenceMatcher(None, s, t, autojunk=autojunk).get_matching_blocks():
            if blk.size >= minblk:
                kc[blk.a:blk.a + blk.size] = True; m[blk.a:blk.a + blk.size] = np.arange(blk.b, blk.b + blk.size)
        return kc, m


    def nw_shared(s, t, minrun=1):
        a = al.align(s, t)[0]; kc = np.zeros(len(s), bool); m = np.full(len(s), -1)
        for (s0, s1), (t0, t1) in zip(*a.aligned):
            for i, j in zip(range(s0, s1), range(t0, t1)):
                if s[i] == t[j]: kc[i] = True; m[i] = j
        if minrun > 1:
            for st, en in runs(kc):
                if en - st + 1 < minrun: kc[st:en + 1] = False; m[st:en + 1] = -1
        return kc, m


    Fc = pd.read_csv(os.path.join(HERE, '..', 'idr_architecture', 'outputs', 'canonical_idr_fate.csv'))
    Fc = Fc[Fc.canonical.isin(B) & Fc.alternative.isin(B)].copy()
    METHODS = ['difflib', 'difflib_autojunk', 'difflib_min5', 'difflib_min10', 'needle', 'needle_min5', 'codon_first_base', 'codon_full', 'codon_full_same_aa']
    rows, match_stats, t0 = [], {m: [0, 0] for m in ['difflib', 'difflib_autojunk', 'difflib_min5', 'needle', 'needle_min5']}, time.time()
    n_skip_nw = 0
    for (c, a), d in Fc.groupby(['canonical', 'alternative']):
        s, t = seq[c], seq[a]
        sh = {'difflib': difflib_shared(s, t, False), 'difflib_autojunk': difflib_shared(s, t, True),
              'difflib_min5': difflib_shared(s, t, False, 5), 'difflib_min10': difflib_shared(s, t, False, 10)}
        if len(s) <= MAXLEN_NW and len(t) <= MAXLEN_NW:
            sh['needle'] = nw_shared(s, t); sh['needle_min5'] = nw_shared(s, t, 5)
        else:
            n_skip_nw += 1
        # do aligned residue pairs share a genomic codon?
        fa = dict(zip(full[a].tolist(), range(len(t))))
        for mth in match_stats:
            if mth not in sh: continue
            kc, m = sh[mth]; ii = np.where(kc)[0]
            match_stats[mth][0] += len(ii); match_stats[mth][1] += int(sum(fa.get(k) != j for k, j in zip(full[c][ii].tolist(), m[ii].tolist())))
        pos = np.array([fa.get(k, -1) for k in full[c].tolist()])
        cod = {'codon_first_base': np.isin(first[c], first[a]), 'codon_full': pos >= 0,
               'codon_full_same_aa': (pos >= 0) & (aa[c] == aa[a][np.maximum(pos, 0)])}
        for r in d.itertuples():
            out = dict(gene=r.gene, canonical=c, alternative=a, idr_len=r.idr_len, original_alignment_fate=r.fate)
            for mth in METHODS:
                x = sh[mth][0] if mth in sh else cod.get(mth)
                out[mth] = classify(x[r.idr_start - 1:r.idr_end]) if x is not None else None
            rows.append(out)
    M = pd.DataFrame(rows); M.to_csv(os.path.join(OUT, 'fate_robustness_methods.csv'), index=False)
    K['R1_runtime_s'] = round(time.time() - t0); K['R1_pairs_skipped_for_needle_length'] = n_skip_nw
    K['R1_original_difflib_reproduced'] = float((M.difflib == M.original_alignment_fate).mean())
    Mc = M.dropna(subset=['needle'])
    K['R1_idrs'] = dict(all=len(M), with_needle=len(Mc), genes=int(Mc.gene.nunique()))
    fates = ['intact', 'truncated', 'removed', 'internal deletion']
    K['R1_fate_pct'] = {m: (100 * Mc[m].value_counts(normalize=True)).reindex(fates).fillna(0).round(2).to_dict() for m in METHODS}
    gw = lambda col: (100 * pd.crosstab(Mc.gene, Mc[col], normalize='index').reindex(columns=fates).fillna(0).mean()).round(2).to_dict()
    K['R1_fate_pct_gene_weighted'] = {m: gw(m) for m in METHODS}
    K['R1_agreement_with_codon_full'] = {m: round(float((Mc[m] == Mc.codon_full).mean()), 4) for m in METHODS}
    K['R1_crosstab_vs_codon_full'] = {m: pd.crosstab(Mc[m], Mc.codon_full).to_dict(orient='index') for m in ['difflib', 'needle', 'needle_min5']}
    K['R1_aligned_identical_residues_not_same_codon'] = {m: dict(aligned=v[0], different_codon=v[1], frac=round(v[1] / v[0], 4) if v[0] else None) for m, v in match_stats.items()}
    # gene bootstrap of the internal-deletion difference (alignment minus codon_full)
    gidx = Mc.groupby('gene').indices; gl = np.array(list(gidx))
    def idd(mth, sample):
        z = Mc.iloc[np.concatenate([gidx[g] for g in sample])]; return 100 * ((z[mth] == 'internal deletion').mean() - (z.codon_full == 'internal deletion').mean())
    K['R1_internal_deletion_excess_vs_codon_full_pp'] = {}
    for mth in ['difflib', 'difflib_autojunk', 'difflib_min5', 'needle', 'needle_min5']:
        bs = [idd(mth, rng.choice(gl, len(gl))) for _ in range(500)]
        K['R1_internal_deletion_excess_vs_codon_full_pp'][mth] = [round(idd(mth, gl), 2), round(np.percentile(bs, 2.5), 2), round(np.percentile(bs, 97.5), 2)]
    print(json.dumps({k: v for k, v in K.items() if k.startswith('R1') or k.startswith('R2')}, indent=1, default=str), flush=True)


# ---------------- R3 + R4: all ordered pairs, full-codon key, disorder retained
fr, ret = [], []
masks = {a: mask(a) for a in B}
for gene, accs in genes.items():
    accs = list(accs); n = len(accs)
    if n < 2: continue
    for a, b in itertools.permutations(accs, 2):
        ib = dict(zip(full[b].tolist(), range(len(full[b]))))
        pos = np.array([ib.get(k, -1) for k in full[a].tolist()]); inB = pos >= 0
        inB1 = np.isin(first[a], first[b])
        for s, e in segs.get(a, []):
            x = inB[s - 1:e]; f = classify(x)
            dB = masks[b][pos[s - 1:e][x]].mean() if x.any() else np.nan
            fr.append((gene, n, a, b, a == can_of.get(gene), b == can_of.get(gene), e - s + 1, f, classify(inB1[s - 1:e]), dB))
if True:
    R = pd.DataFrame(fr, columns=['gene', 'n_iso', 'a', 'b', 'a_can', 'b_can', 'idr_len', 'fate_full', 'fate', 'frac_disordered_in_b'])
    R.to_csv(os.path.join(OUT, 'fate_robustness_all_pairs.csv.gz'), index=False)
K['R3_key_note'] = "R3 'fate' uses the first-base codon key (as in the paper); R3_full_key_* use all three bases"
K['R3_full_key_pooled'] = (100 * R.fate_full.value_counts(normalize=True)).reindex(fates).fillna(0).round(2).to_dict()
K['R3_first_vs_full_key_transitions'] = {f'{a}->{b}': int(n) for (a, b), n in R[R.fate != R.fate_full].groupby(['fate', 'fate_full']).size().items()}
K['R3_idr_comparisons'] = len(R); K['R3_genes'] = int(R.gene.nunique())
pct = lambda z: (100 * z.fate.value_counts(normalize=True)).reindex(fates).fillna(0).round(2).to_dict()
gidx = R.groupby('gene').indices; gl = np.array(list(gidx))
def gene_w(z):
    return (100 * pd.crosstab(z.gene, z.fate, normalize='index').reindex(columns=fates).fillna(0).mean()).round(2).to_dict()
K['R3_pooled'] = pct(R); K['R3_gene_weighted'] = gene_w(R)
bsp, bsg = [], []
for _ in range(300):
    samp = rng.choice(gl, len(gl)); parts = [gidx[g] for g in samp]
    z = R.iloc[np.concatenate(parts)].reset_index(drop=True).assign(gene=np.repeat(np.arange(len(samp)), [len(p) for p in parts]))  # resampled genes stay distinct
    bsp.append(pct(z)); bsg.append(gene_w(z))
K['R3_pooled_ci'] = {f: [round(np.percentile([b[f] for b in bsp], 2.5), 2), round(np.percentile([b[f] for b in bsp], 97.5), 2)] for f in fates}
K['R3_gene_weighted_ci'] = {f: [round(np.percentile([b[f] for b in bsg], 2.5), 2), round(np.percentile([b[f] for b in bsg], 97.5), 2)] for f in fates}
per_gene = R.groupby('gene').size().sort_values(ascending=False)
K['R3_share_of_comparisons_top_1pct_genes'] = round(float(per_gene.head(int(0.01 * len(per_gene))).sum() / len(R)), 3)
K['R3_share_of_comparisons_top_5pct_genes'] = round(float(per_gene.head(int(0.05 * len(per_gene))).sum() / len(R)), 3)
K['R3_top_genes'] = per_gene.head(10).to_dict()
K['R3_without_top_5pct_genes'] = pct(R[~R.gene.isin(per_gene.head(int(0.05 * len(per_gene))).index)])
R['stratum'] = pd.cut(R.n_iso, [1, 2, 3, 5, 10 ** 4], labels=['2', '3', '4-5', '6+'])
K['R3_by_isoform_count'] = {str(s): dict(genes=int(z.gene.nunique()), comparisons=len(z), **pct(z)) for s, z in R.groupby('stratum', observed=True)}
one = R[R.a_can & ~R.b_can]; pick = one.groupby('gene').b.agg(lambda v: rng.choice(v.unique()))
K['R3_one_random_canonical_to_alternative_pair_per_gene'] = dict(genes=len(pick), **pct(one[one.b == one.gene.map(pick)]))
K['R3_by_direction'] = {'canonical_to_alternative': pct(R[R.a_can & ~R.b_can]), 'alternative_to_canonical': pct(R[~R.a_can & R.b_can]),
                        'alternative_to_alternative': pct(R[~R.a_can & ~R.b_can])}
# R4 disorder of retained codons
ok = R.dropna(subset=['frac_disordered_in_b'])
K['R4_retained_codons_disordered_in_b'] = {f: dict(n=int((ok.fate == f).sum()), mean=round(ok[ok.fate == f].frac_disordered_in_b.mean(), 4),
                                                   frac_ge_0_9=round(float((ok[ok.fate == f].frac_disordered_in_b >= 0.9).mean()), 4),
                                                   frac_eq_1=round(float((ok[ok.fate == f].frac_disordered_in_b == 1).mean()), 4),
                                                   frac_lt_0_5=round(float((ok[ok.fate == f].frac_disordered_in_b < 0.5).mean()), 4))
                                           for f in ['intact', 'truncated', 'internal deletion']}
it = ok[ok.fate == 'intact']
K['R4_intact_by_length'] = {lab: round(float((it[(it.idr_len >= lo) & (it.idr_len <= hi)].frac_disordered_in_b >= 0.9).mean()), 4) for lo, hi, lab in LBINS}
json.dump(K, open(os.path.join(OUT, 'fate_robustness.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps({k: v for k, v in K.items() if k.startswith('R3') or k.startswith('R4')}, indent=1, default=str))
