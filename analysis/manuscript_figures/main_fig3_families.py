"""Main Figure 3: protein families and annotations resolve disorder beyond broad class.

a, the ten most and ten least disordered UniProt families (>= 15 genes) with the variance explained by family, class and
shuffled families; b, median canonical %IDR of genes with each annotation vs all other genes. (Disease inheritance moved to
Extended Data Fig. 3, ed3_disease.py.)
"""
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from style import *

PC = AN / 'protein_classes' / 'outputs'
G = pd.read_csv(PC / 'gene_annotations.csv', index_col=0); K = json.load(open(PC / 'key_numbers.json')); CL = K['classes_by_median_pct_idr']
I = pd.read_csv(PC / 'I_families.csv'); J = pd.read_csv(PC / 'J_features.csv')
fig = plt.figure(figsize=(7.1, 5.6)); gs = fig.add_gridspec(1, 2, wspace=0.85)
# a families
ax = fig.add_subplot(gs[0, 0]); sel = pd.concat([I.head(10), I.tail(10)]); fam = G.dropna(subset=['family'])
y = np.arange(len(sel))[::-1]
for yi, r in zip(y, sel.itertuples()):
    v = fam[fam.family == r.family].pct_idr.values; q = np.percentile(v, [10, 25, 50, 75, 90]); col = ORANGE if yi >= 10 else BLUE
    ax.plot([q[0], q[4]], [yi, yi], color=MUTED, lw=0.7); ax.add_patch(plt.Rectangle((q[1], yi - 0.3), q[3] - q[1], 0.6, color=col, alpha=0.8, lw=0))
    ax.plot([q[2], q[2]], [yi - 0.3, yi + 0.3], color=INK, lw=1)
lab = [r.family.replace(' family', '').replace(' superfamily', ' (super)').replace('Krueppel C2H2-type zinc-finger protein', 'C2H2 zinc finger').replace('Class V-like SAM-binding methyltransferase', 'SAM methyltransferase (class V)').replace('Short-chain dehydrogenases/reductases (SDR)', 'Short-chain dehydrogenase (SDR)')[:40] + f' ({r.genes})' for r in sel.itertuples()]
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=5.6); refline(ax, 9.5, 'y', '\u2191 10 most disordered\n\u2193 10 least disordered', ha='right', va='top')
ax.set_xlabel('canonical %IDR\n(box: IQR; line: 10th–90th percentile)'); ax.set_xlim(0, 100); grid(ax, 'x'); tag(ax, 'a', -0.95)
V = K['I_variance_explained_pct_idr']
ax.text(0.98, 0.02, f'variance explained ({V["genes"]:,} genes)\nfamily {100 * V["family"]:.0f}% · class {100 * V["pclass"]:.0f}%\nshuffled families {100 * V["family_null_mean"]:.0f}%',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=5.6, color=MUTED)
# b annotations
ax = fig.add_subplot(gs[0, 1])
NAMES = {'kw_dna_binding': 'DNA-binding', 'kw_transcription': 'transcription regulation', 'many_phospho': '≥5 phosphosites', 'kw_chromatin_regulator': 'chromatin regulator',
         'kw_mrna_processing': 'mRNA processing', 'kw_zinc_finger': 'zinc finger', 'kw_developmental': 'developmental protein', 'kw_rna_binding': 'RNA-binding',
         'loc_nucleus': 'nucleus', 'coiled_coil': 'coiled coil', 'kw_cell_cycle': 'cell cycle', 'loc_synapse': 'synapse', 'disprot': 'in DisProt',
         'many_interactors': '≥10 interactors', 'condensate': 'condensate (CD-CODE)', 'loc_cytoplasm': 'cytoplasm', 'has_pdb': 'has PDB structure',
         'loc_cell_membrane': 'cell membrane', 'disease': 'disease gene', 'kw_kinase': 'kinase', 'signal_peptide': 'signal peptide',
         'transmembrane': 'transmembrane', 'kw_glycoprotein': 'glycoprotein', 'loc_mitochondrion': 'mitochondrion', 'loc_er': 'ER', 'enzyme': 'enzyme (EC)',
         'kw_transport': 'transport', 'loc_secreted': 'secreted', 'kw_receptor': 'receptor', 'kw_ion_channel': 'ion channel'}
jj = J[J.feature.isin(NAMES)].sort_values('median_pct_idr_in').reset_index(drop=True); y = np.arange(len(jj))
for yi, r in zip(y, jj.itertuples()):
    ax.plot([r.median_pct_idr_out, r.median_pct_idr_in], [yi, yi], color=LIGHT, lw=1.2)
    ax.plot(r.median_pct_idr_out, yi, 'o', ms=3, mfc='white', mec=MUTED, mew=0.7)
    ax.plot(r.median_pct_idr_in, yi, 'o', ms=3.5, color=ORANGE if r.median_pct_idr_in > r.median_pct_idr_out else BLUE)
ax.plot(0.02, 1.012, 'o', ms=3.5, color=ORANGE, transform=ax.transAxes, clip_on=False); ax.text(0.045, 1.012, 'genes with the label', transform=ax.transAxes, va='center', fontsize=5.6)
ax.plot(0.6, 1.012, 'o', ms=3, mfc='white', mec=MUTED, mew=0.7, transform=ax.transAxes, clip_on=False); ax.text(0.625, 1.012, 'all other genes', transform=ax.transAxes, va='center', fontsize=5.6)
ax.set_yticks(y); ax.set_yticklabels([f'{NAMES[f]} ({n:,})' for f, n in zip(jj.feature, jj.genes_in)], fontsize=5.6); ax.set_ylim(-0.7, len(jj) - 0.3)
ax.set_xlabel('% of canonical isoform in IDRs\n(median over genes)'); grid(ax, 'x'); tag(ax, 'b', -0.75)
save(fig, 'main/fig3_families')
