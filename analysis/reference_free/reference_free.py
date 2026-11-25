"""Reference-free versions of the IDR-architecture and composition analyses.

No isoform is used as a reference. Every coding residue of every mapped isoform is keyed by the genomic position of
the first base of its codon (which also fixes the reading frame), so the isoforms of a gene share one coordinate
system. For each gene (>= 2 mapped isoforms with a validated coding structure):

  inclusion(position) = fraction of the gene's mapped isoforms that contain the position
  core                = inclusion 1 (present in every isoform); variable = inclusion < 1
  variable classes    = majority (> 1/2), half (= 1/2), minority (< 1/2)
  internal / terminal = a variable position is internal if core positions lie on both sides of it in transcript order
  disordered position = inside a MetaPredict IDR (>= 20 aa) in at least half of the isoforms that contain it

Analyses
  1. Disorder of variable vs core positions (unique positions per gene; gene-paired tests, gene bootstrap)
  2. Composition of disordered positions by inclusion class vs core disorder (gene-clustered linear models)
  3. IDR fates over all ordered isoform pairs (each IDR of isoform A followed into isoform B), pooled and gene-weighted
  4. Disorder classes from core %IDR (the disorder of the shared backbone) and their agreement with canonical classes;
     exon-level splicing-disorder coupling and IDR fates by core class

Writes outputs/key_numbers.json and CSV tables in outputs/.
"""
import os, json, itertools
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'outputs'); os.makedirs(OUT, exist_ok=True)
A = os.path.join(HERE, '..')
rng = np.random.default_rng(20261026)
AA = {'aro': 'FWY', 'acid': 'DE', 'basic': 'KR', 'pro': 'P', 'st': 'ST', 'qn': 'QN', 'gly': 'G', 'hyd': 'AILMV', 'his': 'H', 'cys': 'C'}
FEATS = list(AA) + ['fcr', 'ncpr']
LBINS = [(20, 50, '20-50'), (51, 100, '51-100'), (101, 200, '101-200'), (201, 10 ** 6, '>200')]

E = pd.read_parquet(os.path.join(ROOT, 'data/genomic/processed/exon_level_table.parquet'),
                    columns=['isoform_accession', 'base_accession', 'is_canonical', 'is_coding_exon', 'strand', 'cds_start', 'cds_end', 'aa_start', 'chr'])
E = E[E.is_coding_exon & E.cds_start.notna()]
iso = pd.read_csv(os.path.join(ROOT, 'data/annotated/IDRisoforms_df.csv'), usecols=['isoform_accession', 'Sequence'])
seq = dict(zip(iso.isoform_accession, iso.Sequence))
idr = pd.read_csv(os.path.join(ROOT, 'data/annotated/idr_df.csv')); idr = idr[idr.idr_aaLen >= 20]
segs = {a: [(int(s), int(e)) for s, e in zip(d.idr_start_aa, d.idr_end_aa)] for a, d in idr.groupby('isoform_accession')}

def mask(a):
    m = np.zeros(len(seq[a]), bool)
    for s, e in segs.get(a, []): m[s - 1:e] = True
    return m

def codon_keys(d):
    d = d.sort_values('aa_start', kind='stable'); strand = d.strand.iloc[0]
    parts = [np.arange(int(r.cds_start), int(r.cds_end) + 1) if strand == '+' else np.arange(int(r.cds_end), int(r.cds_start) - 1, -1)
             for r in d.itertuples()]
    b = np.concatenate(parts)
    return b[0::3], len(b)

def runs(x):
    d = np.diff(np.concatenate([[0], x.astype(int), [0]]))
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0] - 1))

def comp(s):
    n = len(s); r = {k: sum(s.count(c) for c in v) / n for k, v in AA.items()}
    r['fcr'] = r['acid'] + r['basic']; r['ncpr'] = r['basic'] - r['acid']; return r

def lbin(L):
    return next(lab for lo, hi, lab in LBINS if lo <= L <= hi)

gene_rows, chunk_rows, fate_counts, pair_rows, incl_rows = [], [], {}, [], []
n_skipped = 0
for gene, d in E.groupby('base_accession'):
    isos = {}
    for a, da in d.groupby('isoform_accession'):
        if a not in seq: continue
        keys, nb = codon_keys(da)
        if nb != 3 * len(seq[a]): n_skipped += 1; continue
        isos[a] = (keys, mask(a))
    n = len(isos)
    if n < 2: continue
    strand = d.strand.iloc[0]
    # ---- gene position table
    allk = np.concatenate([k for k, _ in isos.values()])
    alld = np.concatenate([m for _, m in isos.values()])
    allaa = np.concatenate([np.frombuffer(seq[a].encode(), dtype='S1') for a in isos])
    P = pd.DataFrame({'key': allk, 'dis': alld, 'aa': allaa}).groupby('key').agg(cnt=('dis', 'size'), dis=('dis', 'mean'), aa=('aa', 'first'))
    P['incl'] = P.cnt / n; P['disordered'] = P.dis >= 0.5
    P = P.sort_index(ascending=(strand == '+'))
    core = (P.cnt == n).to_numpy()
    if core.any():
        idx = np.arange(len(P)); first_core, last_core = idx[core][0], idx[core][-1]
        internal = (idx > first_core) & (idx < last_core)
    else:
        internal = np.zeros(len(P), bool)
    P['cls'] = np.select([core, P.incl > 0.5, P.incl == 0.5], ['core', 'majority', 'half'], 'minority')
    P['where'] = np.where(core, 'core', np.where(internal, 'internal', 'terminal'))
    v = P[~core]; c = P[core]
    gene_rows.append(dict(gene=gene, n_iso=n, n_pos=len(P), n_core=int(core.sum()), n_var=int((~core).sum()),
                          dis_core=c.disordered.mean() if len(c) else np.nan, dis_var=v.disordered.mean() if len(v) else np.nan,
                          dis_var_internal=v[v['where'] == 'internal'].disordered.mean() if (v['where'] == 'internal').any() else np.nan,
                          dis_var_terminal=v[v['where'] == 'terminal'].disordered.mean() if (v['where'] == 'terminal').any() else np.nan,
                          core_pct_idr=100 * c.disordered.mean() if len(c) else np.nan,
                          median_iso_pct_idr=float(np.median([100 * m.mean() for _, m in isos.values()]))))
    # inclusion-binned disorder (unique positions)
    for lab, g in P.groupby('cls'):
        incl_rows.append(dict(gene=gene, n_iso=n, cls=lab, n=len(g), n_dis=int(g.disordered.sum())))
    # ---- composition chunks: disordered positions by class (variable split by inclusion and position)
    Dp = P[P.disordered]
    groups = {'core': Dp[Dp.cls == 'core']}
    for lab in ['majority', 'half', 'minority']:
        groups[f'{lab}|internal'] = Dp[(Dp.cls == lab) & (Dp['where'] == 'internal')]
        groups[f'{lab}|terminal'] = Dp[(Dp.cls == lab) & (Dp['where'] == 'terminal')]
    groups['variable|internal'] = Dp[(Dp.cls != 'core') & (Dp['where'] == 'internal')]
    groups['variable|terminal'] = Dp[(Dp.cls != 'core') & (Dp['where'] == 'terminal')]
    for lab, g in groups.items():
        if len(g) >= 10:
            chunk_rows.append(dict(gene=gene, n_iso=n, kind=lab, n=len(g), **comp(b''.join(g.aa.to_list()).decode())))
    # ---- all ordered pairs: IDR fates and spliced vs shared disorder
    keysets = {a: set(k.tolist()) for a, (k, _) in isos.items()}
    for a, b in itertools.permutations(isos, 2):
        ka, ma = isos[a]; inB = np.fromiter((x in keysets[b] for x in ka), bool, len(ka))
        pair_rows.append((gene, int((~inB).sum()), int(ma[~inB].sum()), int(inB.sum()), int(ma[inB].sum())))
        for s, e in segs.get(a, []):
            r = inB[s - 1:e]; f = r.mean()
            if f == 1: cls = 'intact'
            elif f == 0: cls = 'removed'
            else:
                lost = runs(~r); L = len(r)
                cls = 'truncated' if all(st == 0 or en == L - 1 for st, en in lost) and len(lost) <= 2 else 'internal deletion'
            k = (gene, lbin(e - s + 1), cls); fate_counts[k] = fate_counts.get(k, 0) + 1

G = pd.DataFrame(gene_rows); Ch = pd.DataFrame(chunk_rows); IN = pd.DataFrame(incl_rows)
PR = pd.DataFrame(pair_rows, columns=['gene', 'only_a', 'only_a_dis', 'shared', 'shared_dis'])
F = pd.Series(fate_counts).rename_axis(['gene', 'lbin', 'fate']).rename('n').reset_index()
G.to_csv(os.path.join(OUT, 'gene_core_variable.csv'), index=False); Ch.to_csv(os.path.join(OUT, 'disorder_chunks_by_inclusion.csv'), index=False); F.to_csv(os.path.join(OUT, 'all_pairs_idr_fates.csv'), index=False)
K = {'genes': len(G), 'isoforms': int(G.n_iso.sum()), 'isoforms_skipped_length_mismatch': n_skipped,
     'genes_ge3_isoforms': int((G.n_iso >= 3).sum()), 'ordered_pairs': len(PR),
     'positions': int(G.n_pos.sum()), 'core_positions': int(G.n_core.sum()), 'variable_positions': int(G.n_var.sum())}

# 1. variable vs core disorder
pool = IN.groupby('cls')[['n', 'n_dis']].sum(); pool['frac'] = pool.n_dis / pool.n
K['pooled_frac_disordered_by_inclusion'] = pool.frac.round(4).to_dict()
K['pooled_positions_by_inclusion'] = pool.n.astype(int).to_dict()
var = IN[IN.cls != 'core'].groupby('gene')[['n', 'n_dis']].sum(); cor = IN[IN.cls == 'core'].set_index('gene')[['n', 'n_dis']]
K['pooled_frac_disordered'] = dict(core=round(cor.n_dis.sum() / cor.n.sum(), 4), variable=round(var.n_dis.sum() / var.n.sum(), 4))
pg = G.dropna(subset=['dis_core', 'dis_var'])
K['gene_paired'] = dict(genes=len(pg), mean_core=round(pg.dis_core.mean(), 4), mean_variable=round(pg.dis_var.mean(), 4),
                        frac_var_higher=round(float((pg.dis_var > pg.dis_core).mean()), 3), frac_var_lower=round(float((pg.dis_var < pg.dis_core).mean()), 3),
                        wilcoxon_p=float(f'{stats.wilcoxon(pg.dis_var, pg.dis_core).pvalue:.3g}'))
idx = {g: i for i, g in enumerate(IN.gene.unique())}
gv = IN.assign(var=IN.cls != 'core').groupby(['gene', 'var'])[['n', 'n_dis']].sum().unstack().dropna()
def boot_diff(sample):
    s = gv.loc[sample]; return s[('n_dis', True)].sum() / s[('n', True)].sum() - s[('n_dis', False)].sum() / s[('n', False)].sum()
genes = gv.index.to_numpy()
bs = [boot_diff(rng.choice(genes, len(genes))) for _ in range(1000)]
K['pooled_diff_variable_minus_core'] = [round(boot_diff(genes), 4), round(np.percentile(bs, 2.5), 4), round(np.percentile(bs, 97.5), 4)]
for w in ['internal', 'terminal']:
    q = G.dropna(subset=['dis_core', f'dis_var_{w}'])
    K[f'gene_paired_{w}'] = dict(genes=len(q), mean_variable=round(q[f'dis_var_{w}'].mean(), 4), mean_core=round(q.dis_core.mean(), 4),
                                 wilcoxon_p=float(f'{stats.wilcoxon(q[f"dis_var_{w}"], q.dis_core).pvalue:.3g}'))
# genes with >= 3 isoforms: minority vs majority inclusion
g3 = IN[IN.n_iso >= 3].pivot_table(index='gene', columns='cls', values=['n', 'n_dis'], aggfunc='sum')
m3 = g3.dropna(subset=[('n', 'minority'), ('n', 'majority')])
fr = lambda c: m3[('n_dis', c)] / m3[('n', c)]
K['ge3_minority_vs_majority'] = dict(genes=len(m3), mean_minority=round(fr('minority').mean(), 4), mean_majority=round(fr('majority').mean(), 4),
                                     mean_core=round(fr('core').mean(), 4) if ('n', 'core') in m3 else None,
                                     wilcoxon_p=float(f'{stats.wilcoxon(fr("minority"), fr("majority")).pvalue:.3g}'))

# 2. composition vs core disorder
crow = []
for kind in [k for k in Ch.kind.unique() if k != 'core']:
    x = Ch[Ch.kind.isin([kind, 'core'])].copy(); x['s'] = (x.kind == kind).astype(int); x['logn'] = np.log10(x.n)
    if x.s.sum() < 50: continue
    for f in FEATS:
        m = smf.ols(f'{f} ~ s + logn', x).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(x.gene)[0]})
        crow.append(dict(kind=kind, feature=f, n_chunks=int(x.s.sum()), diff_pp=100 * m.params.s, lo=100 * m.conf_int().loc['s', 0],
                         hi=100 * m.conf_int().loc['s', 1], p=m.pvalues.s))
C = pd.DataFrame(crow); C['q_BH'] = C.groupby('kind').p.transform(lambda p: multipletests(p, method='fdr_bh')[1])
C.to_csv(os.path.join(OUT, 'composition_vs_core.csv'), index=False)
K['composition_vs_core'] = {k: {r.feature: [round(r.diff_pp, 2), round(r.lo, 2), round(r.hi, 2), float(f'{r.q_BH:.2g}')] for r in g.itertuples()} | {'n_chunks': int(g.n_chunks.iloc[0])}
                            for k, g in C.groupby('kind')}

# 3. all-pairs IDR fates
fates = ['intact', 'internal deletion', 'truncated', 'removed']
tot = F.groupby('fate').n.sum().reindex(fates).fillna(0); K['all_pairs_fate_pct_pooled'] = (100 * tot / tot.sum()).round(1).to_dict()
gw = F.pivot_table(index='gene', columns='fate', values='n', aggfunc='sum').reindex(columns=fates).fillna(0)
K['all_pairs_fate_pct_gene_weighted'] = (100 * gw.div(gw.sum(axis=1), axis=0).mean()).round(1).to_dict()
K['all_pairs_idr_comparisons'] = int(tot.sum())
lb = F.pivot_table(index='lbin', columns='fate', values='n', aggfunc='sum').reindex(columns=fates).fillna(0)
K['all_pairs_fate_pct_by_length'] = (100 * lb.div(lb.sum(axis=1), axis=0)).round(1).to_dict(orient='index')
PRg = PR.groupby('gene').sum()
K['all_pairs_disorder'] = dict(pooled_isoform_specific=round(PR.only_a_dis.sum() / PR.only_a.sum(), 4), pooled_shared=round(PR.shared_dis.sum() / PR.shared.sum(), 4))

# 4. core-%IDR disorder classes
G['core_class'] = np.select([G.core_pct_idr >= 75, G.core_pct_idr <= 15], ['high', 'very low'], 'intermediate')
G.loc[G.core_pct_idr.isna(), 'core_class'] = np.nan
can = pd.read_csv(os.path.join(A, 'idr_landscape', 'outputs', 'gene_disorder_class.csv'), index_col=0).dclass
G['canonical_class'] = G.gene.map(can)
ag = G.dropna(subset=['core_class', 'canonical_class'])
K['class_agreement'] = dict(genes=len(ag), agree=round(float((ag.core_class == ag.canonical_class).mean()), 3),
                            crosstab=pd.crosstab(ag.canonical_class, ag.core_class).to_dict(orient='index'),
                            core_class_sizes=G.core_class.value_counts().to_dict(), genes_without_core=int(G.core_pct_idr.isna().sum()))
ux = pd.read_csv(os.path.join(A, 'reanalysis_2026-10-26', 'outputs', 'unique_exons.csv'))
ux['core_class'] = ux.gene.map(G.set_index('gene').core_class)
uidx = ux.groupby('gene').indices
def ediff(rows):
    z = ux.iloc[rows]; return z[z.alt == 1].idr_frac.mean() - z[z.alt == 0].idr_frac.mean()
K['exon_alt_minus_const_by_core_class'] = {}
for cl in ['very low', 'intermediate', 'high']:
    gs = ux[ux.core_class == cl].gene.unique()
    if len(gs) < 20: continue
    b = [ediff(np.concatenate([uidx[g] for g in rng.choice(gs, len(gs))])) for _ in range(1000)]
    K['exon_alt_minus_const_by_core_class'][cl] = dict(genes=int(len(gs)), diff=round(ediff(np.concatenate([uidx[g] for g in gs])), 4),
                                                        lo=round(np.percentile(b, 2.5), 4), hi=round(np.percentile(b, 97.5), 4))
# classes from the median %IDR over a gene's isoforms (also reference-free). Core %IDR is not used for this test:
# it is built from the same positions as the constitutive exons, so selecting genes on it biases alt - const.
G['median_class'] = np.select([G.median_iso_pct_idr >= 75, G.median_iso_pct_idr <= 15], ['high', 'very low'], 'intermediate')
ux['median_class'] = ux.gene.map(G.set_index('gene').median_class)
K['exon_alt_minus_const_by_median_class'] = {}
for cl in ['very low', 'intermediate', 'high']:
    gs = ux[ux.median_class == cl].gene.unique()
    b = [ediff(np.concatenate([uidx[g] for g in rng.choice(gs, len(gs))])) for _ in range(1000)]
    K['exon_alt_minus_const_by_median_class'][cl] = dict(genes=int(len(gs)), diff=round(ediff(np.concatenate([uidx[g] for g in gs])), 4),
                                                          lo=round(np.percentile(b, 2.5), 4), hi=round(np.percentile(b, 97.5), 4))
K['median_class_sizes'] = G.median_class.value_counts().to_dict()
F['core_class'] = F.gene.map(G.set_index('gene').core_class)
fc = F.dropna(subset=['core_class']).pivot_table(index='core_class', columns='fate', values='n', aggfunc='sum').reindex(columns=fates).fillna(0)
K['all_pairs_fate_pct_by_core_class'] = (100 * fc.div(fc.sum(axis=1), axis=0)).round(1).to_dict(orient='index')

# canonical-referenced numbers for comparison
arch = json.load(open(os.path.join(A, 'idr_architecture', 'outputs', 'key_numbers.json')))
K['canonical_referenced_for_comparison'] = dict(fate_pct=arch.get('canonical_idr_fate_pct'), pooled_disorder=arch.get('pooled_fraction_disordered'))
json.dump(K, open(os.path.join(OUT, 'key_numbers.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps({k: v for k, v in K.items() if k != 'composition_vs_core'}, indent=1, default=str))
print(C[C.q_BH < 0.05].round(3).to_string())
