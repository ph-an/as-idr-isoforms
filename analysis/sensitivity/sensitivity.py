"""Sensitivity of the main results to how disorder is called.

Disorder calls compared (each gives a per-residue disordered/not mask for every cohort isoform):
  stored        IDRs used throughout the paper: MetaPredict v3.1 predict_disorder_domains (threshold 0.5, smoothing, gap
                closure, minimum IDR 12) and then segments >= 20 residues (data/annotated/idr_df.csv)
  v3_dom_t0.4 / v3_dom_t0.6     the same domain call at disorder thresholds 0.4 and 0.6, segments >= 20
  v3_dom_min10 / min30 / min50  threshold 0.5, minimum segment length 10, 30 or 50
  v3_raw_t0.5   raw per-residue scores >= 0.5 with no smoothing, runs >= 20 (the rule the Methods text described)
The v3 calls need per-residue scores from predict_scores.py; if those are missing only stored is run.

For each call we recompute
  S1  median canonical %IDR per protein class (and the Spearman correlation of the class ordering with the stored call)
  S2  alternative minus constitutive mean exon IDR fraction: raw (gene bootstrap) and adjusted (log length, exon position,
      log transcripts; gene-clustered errors); number of classes whose bootstrap interval excludes zero
  S3  within-gene effect of a lineage-restricted exon (splice site conserved in < 100 of 405 genomes; gene fixed effects)
  S4  change in canonical %IDR per 100 Myr of gene age among family genes, without and with family fixed effects
Writes outputs/sensitivity.json and outputs/sensitivity_summary.csv.
"""
import os, glob, pickle, json, sys
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); A = os.path.join(ROOT, 'analysis')
OUT = os.path.join(HERE, 'outputs'); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261105)
iso = pd.read_csv(os.path.join(ROOT, 'data/annotated/IDRisoforms_df.csv'), usecols=['isoform_accession', 'base_accession', 'is_canonical', 'Sequence'])
seq = dict(zip(iso.isoform_accession, iso.Sequence))
idr = pd.read_csv(os.path.join(ROOT, 'data/annotated/idr_df.csv')); idr = idr[idr.idr_aaLen >= 20]


def runs(x):
    d = np.diff(np.concatenate([[0], x.astype(int), [0]]))
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))


def keep_runs(x, L):
    m = np.zeros(len(x), bool)
    for s, e in runs(x):
        if e - s >= L: m[s:e] = True
    return m


MASKS = {}
st = {a: np.zeros(len(s), bool) for a, s in seq.items()}
for r in idr.itertuples():
    if r.isoform_accession in st: st[r.isoform_accession][r.idr_start_aa - 1:r.idr_end_aa] = True
MASKS['stored'] = st
# MetaPredict v3 variants from stored per-residue scores
parts = sorted(glob.glob(os.path.join(HERE, 'outputs', 'scores', 'part_*.pkl')))
if parts:
    from metapredict.backend.domain_definition import get_domains
    SC = {}
    for f in parts: SC.update(pickle.load(open(f, 'rb')))
    if len(SC) < len(seq) * 0.99: print(f'only {len(SC)} of {len(seq)} isoforms have scores; v3 variants skipped'); SC = {}
    if SC:
        def dom(a, t, L):
            m = np.zeros(len(seq[a]), bool)
            for s, e, *_ in get_domains(seq[a], SC[a].astype(float), disorder_threshold=t)[1]:
                if e - s >= L: m[s:e] = True
            return m
        for name, t, L in [('v3_dom_t0.5_check', 0.5, 20), ('v3_dom_t0.4', 0.4, 20), ('v3_dom_t0.6', 0.6, 20), ('v3_dom_min10', 0.5, 10),
                           ('v3_dom_min30', 0.5, 30), ('v3_dom_min50', 0.5, 50)]:
            MASKS[name] = {a: dom(a, t, L) for a in seq}
        MASKS['v3_raw_t0.5'] = {a: keep_runs(SC[a] >= 0.5, 20) for a in seq}
        chk = MASKS['v3_dom_t0.5_check']
        print('residue agreement of recomputed default call with stored IDRs:', np.mean(np.concatenate([chk[a] == st[a] for a in seq])), flush=True)

# MetaPredict V1 (trained on MobiDB consensus, not on pLDDT) for the subset run by predict_v1_subset.py
v1parts = sorted(glob.glob(os.path.join(HERE, 'outputs', 'v1', 'part_*.pkl')) + glob.glob(os.path.join(HERE, 'outputs', 'v1_all', 'part_*.pkl')))
if v1parts:
    from metapredict.backend.domain_definition import get_domains
    V1 = {}
    for f in v1parts: V1.update(pickle.load(open(f, 'rb')))
    V1 = {a: v for a, v in V1.items() if a in seq and len(v) == len(seq[a])}
    def dom1(a, t=0.3, L=20):
        m = np.zeros(len(seq[a]), bool)
        for s_, e_, *_ in get_domains(seq[a], V1[a].astype(float), disorder_threshold=t)[1]:
            if e_ - s_ >= L: m[s_:e_] = True
        return m
    MASKS['v1_dom_t0.3'] = {a: dom1(a) for a in V1}
    MASKS['v1_raw_t0.3'] = {a: keep_runs(V1[a] >= 0.3, 20) for a in V1}
    MASKS['stored|v1_subset'] = {a: st[a] for a in V1}
# MobiDB-lite consensus ('Disordered' regions in UniProt), canonical isoforms only
U = pd.read_csv(os.path.join(ROOT, 'data/raw/uniprot/uniprotkb_annotated'), sep='\t', usecols=['Entry', 'Region', 'Sequence'], low_memory=False).dropna(subset=['Sequence'])
import re
mob = {}
can_set = set(iso[iso.is_canonical].isoform_accession)
for e, reg, sq_ in zip(U.Entry, U.Region.fillna(''), U.Sequence):
    a = e if e in seq else (e + '-1' if e + '-1' in seq else None)
    if a is None or a not in can_set or seq[a] != sq_: continue
    m = np.zeros(len(sq_), bool)
    for x in re.finditer(r'REGION (\d+)\.\.(\d+); /note="Disordered"', reg): m[int(x.group(1)) - 1:int(x.group(2))] = True
    mob[a] = m
MASKS['mobidb_lite|canonical'] = mob
MASKS['stored|canonical'] = {a: st[a] for a in mob}
print('mobidb canonical isoforms', len(mob), flush=True)

# ---- shared tables
G = pd.read_csv(os.path.join(A, 'protein_classes/outputs/gene_annotations.csv'), index_col=0)
can = iso[iso.is_canonical].drop_duplicates('base_accession').set_index('base_accession').isoform_accession
EX = pd.read_csv(os.path.join(ROOT, 'data/genomic/analysis/exon_idr_annotation.csv'),
                 usecols=['gene', 'iso', 'exon_id', 'exon_rank', 'n_coding_exons', 'aa_start', 'aa_end', 'usage', 'n_tx_gene'])
EX = EX[EX.usage.isin(['alternative', 'constitutive']) & EX.iso.isin(seq.keys())].copy()
CONS = pd.read_csv(os.path.join(A, 'exon_age/outputs/exon_conservation.csv'), usecols=['gene', 'exon_id', 'cons', 'n_sites', 'pos'])
FAMG = G.dropna(subset=['family', 'gene_age']); FAMG = FAMG[FAMG.groupby('family').family.transform('size') >= 15]


def exon_table(M):
    f = [M[i][s - 1:e].mean() if (i in M and e >= s) else np.nan for i, s, e in zip(EX.iso, EX.aa_start.astype(int), EX.aa_end.astype(int))]
    z = EX.assign(idr_frac=f).groupby(['gene', 'exon_id']).agg(idr_frac=('idr_frac', 'mean'), usage=('usage', 'first'), aa_len=('aa_end', 'median'),
                                                              a0=('aa_start', 'median'), first=('exon_rank', 'min'), n_tx=('n_tx_gene', 'first')).reset_index()
    z['aa_len'] = (z.aa_len - z.a0 + 1).clip(lower=1)
    z['alt'] = (z.usage == 'alternative').astype(int); z['log_len'] = np.log10(z.aa_len); z['log_ntx'] = np.log2(z.n_tx)
    z = z.merge(CONS, on=['gene', 'exon_id'], how='left'); z['pclass'] = z.gene.map(G.pclass)
    return z.dropna(subset=['idr_frac'])


def clus(f, d):
    return smf.ols(f, d).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(d.gene)[0]})


def fe(df, y, xs):
    z = df.dropna(subset=[y] + xs).copy(); z = z[z.groupby('gene')[y].transform('size') > 1]
    for c in [y] + xs: z[c + '_d'] = z[c] - z.groupby('gene')[c].transform('mean')
    return smf.ols(f'{y}_d ~ ' + ' + '.join(c + '_d' for c in xs) + ' - 1', z).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(z.gene)[0]})


def boot_alt(z, B=200):
    idx = z.groupby('gene').indices; gl = np.array(list(idx))
    def d(s):
        q = z.iloc[np.concatenate([idx[g] for g in s])]; return q[q.alt == 1].idr_frac.mean() - q[q.alt == 0].idr_frac.mean()
    b = [d(rng.choice(gl, len(gl))) for _ in range(B)]; return d(gl), np.percentile(b, 2.5), np.percentile(b, 97.5)


K, base_rank = {}, None
ONLY = sys.argv[1:]
for name, M in MASKS.items():
    if ONLY and name not in ONLY and name != 'stored': continue
    r = {}
    pct = pd.Series({g: 100 * M[a].mean() for g, a in can.items() if a in M}); G['p'] = pct
    med = G.groupby('pclass').p.median().sort_values(); r['S1_class_medians'] = med.round(1).to_dict()
    r['S1_genome_median_pct'] = round(float(pct.median()), 2)
    if base_rank is None: base_rank = med
    r['S1_spearman_with_stored'] = round(stats.spearmanr(base_rank, med.reindex(base_rank.index)).statistic, 3)
    z = exon_table(M)
    d, lo, hi = boot_alt(z); r['S2_alt_minus_const'] = [round(d, 4), round(lo, 4), round(hi, 4)]
    m = clus('idr_frac ~ alt + log_len + C(np.where(first == 1, "first", "other")) + log_ntx', z)
    r['S2_adjusted'] = [round(m.params.alt, 4), round(m.conf_int().loc['alt', 0], 4), round(m.conf_int().loc['alt', 1], 4)]
    npos = 0
    for c, zc in z.groupby('pclass'):
        if zc.gene.nunique() < 30: continue
        _, l, _ = boot_alt(zc, 100); npos += l > 0
    r['S2_classes_ci_above_0'] = int(npos)
    I = z[(z.pos == 'internal') & (z.n_sites == 2)].copy(); I['young'] = (I.cons < 100).astype(float)
    m = fe(I, 'idr_frac', ['young', 'log_len', 'alt', 'log_ntx'])
    r['S3_young_exon_within_gene'] = [round(m.params.young_d, 4), round(m.conf_int().loc['young_d', 0], 4), round(m.conf_int().loc['young_d', 1], 4), float(f'{m.pvalues.young_d:.2g}')]
    F = FAMG.assign(p=pct.reindex(FAMG.index), age100=FAMG.gene_age / 100).dropna(subset=['p'])
    m0 = smf.ols('p ~ age100', F).fit(cov_type='HC1'); m1 = smf.ols('p ~ age100 + C(family)', F).fit(cov_type='HC1')
    r['S4_slope_per_100Myr'] = [round(m0.params.age100, 3), round(m0.conf_int().loc['age100', 0], 3), round(m0.conf_int().loc['age100', 1], 3)]
    r['S4_slope_family_FE'] = [round(m1.params.age100, 3), round(m1.conf_int().loc['age100', 0], 3), round(m1.conf_int().loc['age100', 1], 3), float(f'{m1.pvalues.age100:.2g}')]
    ga = G.dropna(subset=['gene_age', 'p'])
    r['S5_spearman_age_all_genes'] = [round(v, 4) for v in stats.spearmanr(ga.p, ga.gene_age)] + [len(ga)]
    r['S5_median_by_age_bin'] = ga.groupby('age_bin').p.median().round(1).to_dict()
    r['S5_spearman_by_class'] = {c: round(stats.spearmanr(d.p, d.gene_age).statistic, 3) for c, d in ga.groupby('pclass')}
    K[name] = r; print(name, json.dumps(r), flush=True)
json.dump(K, open(os.path.join(OUT, f'sensitivity{"_" + "_".join(ONLY).replace("|", "-") if ONLY else ""}.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
rows = [dict(call=n, genome_median_pct=r['S1_genome_median_pct'], class_rank_rho=r['S1_spearman_with_stored'], alt_minus_const=r['S2_alt_minus_const'][0],
             alt_lo=r['S2_alt_minus_const'][1], alt_hi=r['S2_alt_minus_const'][2], alt_adjusted=r['S2_adjusted'][0], classes_positive=r['S2_classes_ci_above_0'],
             young_exon=r['S3_young_exon_within_gene'][0], young_lo=r['S3_young_exon_within_gene'][1], young_hi=r['S3_young_exon_within_gene'][2],
             age_rho_all=r['S5_spearman_age_all_genes'][0], age_slope=r['S4_slope_per_100Myr'][0], age_slope_family_fe=r['S4_slope_family_FE'][0], fe_lo=r['S4_slope_family_FE'][1], fe_hi=r['S4_slope_family_FE'][2])
        for n, r in K.items()]
pd.DataFrame(rows).to_csv(os.path.join(OUT, f'sensitivity_summary{"_" + "_".join(ONLY).replace("|", "-") if ONLY else ""}.csv'), index=False); print(pd.DataFrame(rows).to_string())
