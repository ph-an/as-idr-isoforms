"""Example genes: candidate search and side-by-side gene/protein diagrams.

  python ed5_examples.py candidates          # writes outputs/example_candidates.csv
  python ed5_examples.py plot GENE1 GENE2    # writes manuscript/figures/extended_data/ed5_examples.{pdf,png}

Each isoform is drawn as its coding exons in genomic order with introns compressed to a fixed gap, exon width
proportional to coding length, and the residues inside MetaPredict IDRs (>= 20 aa) shaded. Exons are keyed by
coding coordinates, and coding exons present in some but not all of the drawn isoforms are outlined.
"""
import os, sys, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent; ROOT = Path(os.environ.get('TFIDR_ROOT', HERE.parents[1]))
FIG = HERE.parents[1] / 'manuscript' / 'figures'; OUT = HERE / 'outputs'; OUT.mkdir(exist_ok=True)
E = pd.read_parquet(ROOT / 'data/genomic/processed/exon_level_table.parquet')
E = E[E.is_coding_exon & E.aa_start.notna()]
idr = pd.read_csv(ROOT / 'data/annotated/idr_df.csv'); idr = idr[idr.idr_aaLen >= 20]
segs = {a: list(zip(d.idr_start_aa, d.idr_end_aa)) for a, d in idr.groupby('isoform_accession')}
names = pd.read_csv(ROOT / 'data/annotated/IDRisoforms_df.csv', usecols=['base_accession', 'Protein names']).drop_duplicates('base_accession').set_index('base_accession')['Protein names']

def candidates():
    rows = []
    for gid, d in E.groupby('gene_id'):
        isos = d.isoform_accession.unique()
        if not (2 <= len(isos) <= 4) or not d.is_canonical.any(): continue
        can = d[d.is_canonical]
        L = can.aa_end.max()
        if not (250 <= L <= 1000): continue
        u = d.drop_duplicates('exon_id')
        alt = u[(u.exon_usage == 'alternative') & (u.aa_length >= 25)]
        if alt.empty or len(u) > 20: continue
        best = alt.sort_values('idr_frac', ascending=False).iloc[0]
        ordered_can = (can.aa_length * (1 - can.idr_frac)).sum()
        rows.append(dict(gene_id=gid, gene=d.gene_name.iloc[0], acc=can.base_accession.iloc[0], n_iso=len(isos), length=int(L),
                         n_exons=len(u), best_alt_exon_idr=round(best.idr_frac, 2), best_alt_exon_aa=int(best.aa_length),
                         best_alt_first_last=bool(best.is_first_exon or best.is_last_exon), canonical_ordered_aa=int(ordered_can),
                         canonical_idr_frac=round((can.aa_length * can.idr_frac).sum() / can.aa_length.sum(), 2),
                         protein=str(names.get(can.base_accession.iloc[0], ''))[:60]))
    C = pd.DataFrame(rows)
    C = C[(C.best_alt_exon_idr >= 0.9) & (C.canonical_ordered_aa >= 120) & C.canonical_idr_frac.between(0.2, 0.7)]
    C.sort_values(['best_alt_first_last', 'best_alt_exon_aa'], ascending=[True, False]).to_csv(OUT / 'example_candidates.csv', index=False)
    print(len(C), 'candidates ->', OUT / 'example_candidates.csv')

def draw(ax, gene):
    d = E[E.gene_name == gene].copy()
    # key exons by coding coordinates: identical coding segments with different UTRs carry different Ensembl exon IDs
    d['exon_id'] = d.cds_start.astype(int).astype(str) + '-' + d.cds_end.astype(int).astype(str)
    n_iso = d.isoform_accession.nunique()
    d['exon_usage'] = np.where(d.groupby('exon_id').isoform_accession.transform('nunique') < n_iso, 'alternative', 'constitutive')
    ex = d.drop_duplicates('exon_id').sort_values('cds_start').copy(); ex['start'] = ex.cds_start
    strand = ex.strand.iloc[0]
    order = list(ex.exon_id) if strand == '+' else list(ex.exon_id)[::-1]
    width = {e: max(int(r.aa_length), 3) for e, r in ex.set_index('exon_id').iterrows()}
    GAP = max(8, int(0.03 * sum(width.values()))); x0 = {}; x = 0
    for e in order: x0[e] = x; x += width[e] + GAP
    isos = d.groupby('isoform_accession').is_canonical.first().sort_values(ascending=False).index.tolist()
    for i, iso in enumerate(isos):
        y = -i * 1.0; di = d[d.isoform_accession == iso].set_index('exon_id')
        es = [e for e in order if e in di.index]
        ax.plot([x0[es[0]], x0[es[-1]] + width[es[-1]]], [y, y], color='#b9b8b3', lw=0.8, zorder=0)
        sg = segs.get(iso, [])
        for e in es:
            r = di.loc[e]; a0, a1 = int(r.aa_start), int(r.aa_end); w = width[e]
            alt = r.exon_usage == 'alternative'
            ax.add_patch(Rectangle((x0[e], y - 0.32), w, 0.64, fc='#e6e5e0', ec='#eb6834' if alt else 'none', lw=1.4 if alt else 0, zorder=2))
            for s0, s1 in sg:   # IDR residues inside this exon, placed in protein order within the exon box
                lo, hi = max(s0, a0), min(s1, a1)
                if hi < lo: continue
                f0 = (lo - a0) / max(a1 - a0 + 1, 1); f1 = (hi - a0 + 1) / max(a1 - a0 + 1, 1)
                if strand == '-': pass  # exon boxes are drawn in transcript order already (order reversed for '-')
                ax.add_patch(Rectangle((x0[e] + f0 * w, y - 0.32), (f1 - f0) * w, 0.64, fc='#2a78d6', ec='none', zorder=3))
        lab = iso + (' (canonical)' if d[d.isoform_accession == iso].is_canonical.iloc[0] else '')
        pct = 100 * sum(e1 - e0 + 1 for e0, e1 in sg) / int(d[d.isoform_accession == iso].aa_end.max())
        ax.text(-GAP, y, f'{lab}\n{int(d[d.isoform_accession == iso].aa_end.max())} aa, {pct:.0f}% IDR', ha='right', va='center', fontsize=6.3, color='#52514e')
    ax.set_xlim(-GAP * 0.5, x); ax.set_ylim(-len(isos) + 0.4, 0.8); ax.axis('off')
    ax.text(0, 0.62, f'{gene}  ({ex.chr.iloc[0]}, {strand} strand)', fontsize=8, fontweight='bold', va='bottom')

if __name__ == '__main__':
    if sys.argv[1] == 'candidates': candidates()
    else:
        genes = sys.argv[2:4]
        plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'pdf.fonttype': 42})
        fig, axs = plt.subplots(1, 2, figsize=(7.1, 1.9), gridspec_kw={'wspace': 0.55})
        for ax, gname, tag in zip(axs, genes, 'ab'):
            draw(ax, gname); ax.text(-0.02, 1.08, tag, transform=ax.transAxes, fontweight='bold', fontsize=9)
        h = [Rectangle((0, 0), 1, 1, fc='#2a78d6'), Rectangle((0, 0), 1, 1, fc='#e6e5e0'), Rectangle((0, 0), 1, 1, fc='white', ec='#eb6834', lw=1.4)]
        fig.legend(h, ['residues in IDRs', 'ordered residues', 'alternative exon'], loc='lower center', ncol=3, frameon=False, fontsize=6.5, bbox_to_anchor=(0.6, -0.08))
        (FIG / 'extended_data').mkdir(parents=True, exist_ok=True)
        for e in ['pdf', 'png']: fig.savefig(FIG / 'extended_data' / f'ed5_examples.{e}', dpi=300, bbox_inches='tight')
        print('written', FIG / 'extended_data' / 'ed5_examples.pdf')
