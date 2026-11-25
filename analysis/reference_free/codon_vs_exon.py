"""Exon labels versus codon labels (Fig. 2b).

For each gene with at least two mapped isoforms (the same genes as reference_free.py), every unique genomic codon position
is labelled twice: shared or isoform-specific (codon level), and in a constitutive or an alternative exon (exon level; an
exon is alternative if its Ensembl identifier is in some but not all of the gene's mapped transcripts). A codon counts as
in an alternative exon if any isoform that contains it encodes it from an alternative exon. Shared codons inside
alternative exons are split by why the exon is alternative: the coding part is identical in every isoform and only the
untranslated end (and so the exon identifier) differs, or the exon's coding boundary shifts in some isoforms and this
part is kept in all. Isoform-specific codons inside constitutive exons are split into the same DNA read in another frame
or DNA not present in every isoform.

Writes outputs/codon_vs_exon.json.
"""
import os, json, collections
import numpy as np, pandas as pd

ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')
E = pd.read_parquet(os.path.join(ROOT, 'data/genomic/processed/exon_level_table.parquet'),
                    columns=['isoform_accession', 'base_accession', 'is_coding_exon', 'strand', 'cds_start', 'cds_end', 'aa_start', 'exon_usage'])
E = E[E.is_coding_exon & E.cds_start.notna()]
iso = pd.read_csv(os.path.join(ROOT, 'data/annotated/IDRisoforms_df.csv'), usecols=['isoform_accession', 'Sequence'])
L = dict(zip(iso.isoform_accession, iso.Sequence.str.len()))
C = collections.Counter(); W = collections.Counter(); genes = 0
for g, d in E.groupby('base_accession'):
    isos = {}
    for a, da in d.groupby('isoform_accession'):
        if a not in L: continue
        da = da.sort_values('aa_start', kind='stable'); s = da.strand.iloc[0]; b, lab = [], []
        for r in da.itertuples():
            x = np.arange(int(r.cds_start), int(r.cds_end) + 1) if s == '+' else np.arange(int(r.cds_end), int(r.cds_start) - 1, -1)
            b.append(x); lab += [(r.exon_usage == 'alternative', (int(r.cds_start), int(r.cds_end)))] * len(x)
        b = np.concatenate(b)
        if len(b) != 3 * L[a]: continue
        isos[a] = (b[0::3].tolist(), lab[0::3], set(b.tolist()))
    n = len(isos)
    if n < 2: continue
    genes += 1; cnt = collections.Counter(); info = collections.defaultdict(list)
    for k, lab, _ in isos.values():
        for kk, l in zip(k, lab): cnt[kk] += 1; info[kk].append(l)
    shared_bases = set.intersection(*[x[2] for x in isos.values()])
    for kk, c in cnt.items():
        alt = any(l[0] for l in info[kk]); shared = c == n
        C[('shared' if shared else 'isoform_specific') + '|' + ('alternative' if alt else 'constitutive')] += 1
        if shared and alt:
            W['shared_in_alt|same coding part, untranslated end differs' if len({l[1] for l in info[kk]}) == 1 else 'shared_in_alt|coding boundary shifts'] += 1
        if not shared and not alt:
            W['specific_in_const|same DNA, other frame' if kk in shared_bases else 'specific_in_const|DNA not in every isoform'] += 1
spec = C['isoform_specific|alternative'] + C['isoform_specific|constitutive']; inalt = C['isoform_specific|alternative'] + C['shared|alternative']
K = dict(genes=genes, positions=sum(C.values()), counts=dict(C), reasons=dict(W),
         pct_isoform_specific_in_alternative_exons=round(100 * C['isoform_specific|alternative'] / spec, 1),
         pct_alternative_exon_positions_shared=round(100 * C['shared|alternative'] / inalt, 1),
         pct_shared_in_alt_boundary_shift=round(100 * W['shared_in_alt|coding boundary shifts'] / C['shared|alternative'], 1))
json.dump(K, open(os.path.join(OUT, 'codon_vs_exon.json'), 'w'), indent=1); print(json.dumps(K, indent=1))
