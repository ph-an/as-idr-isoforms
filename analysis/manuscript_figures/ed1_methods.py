"""Extended Data Figure 1: methods workflow (redrawn from the 2026 QCB SURP poster and updated).

Eight steps from reviewed UniProt isoforms to the analyses, with the number of isoforms, genes, exons or positions that
pass each step. Step 4 shows how many isoforms are kept at each transcript-matching stage as bars instead of a funnel.
All counts are read from the pipeline outputs.
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, Ellipse, FancyArrowPatch, Polygon
from style import *

GREEN, GREEN_E, PINK, PINK_E = '#bfe0a8', '#6aa84f', '#f6b9b9', '#e06666'
CC = json.load(open(AN / 'manuscript_figures/outputs/cohort_counts.json')); RF = json.load(open(AN / 'reference_free/outputs/key_numbers.json'))
PC = json.load(open(AN / 'protein_classes/outputs/key_numbers.json')); EV = json.load(open(AN / 'evolution/outputs/key_numbers.json'))
EA = json.load(open(AN / 'exon_age/outputs/key_numbers.json'))

fig, ax = plt.subplots(figsize=(7.1, 5.4)); ax.set_xlim(0, 20); ax.set_ylim(0, 14); ax.axis('off')
STEP_COL = ['#3f7fbf', '#7b5ea7', DIS_E, '#2a9d8f', '#c99a06', '#d6619e', '#8c564b', '#3a923a']
W, H = 4.55, 6.3; X = [0.2, 5.2, 10.2, 15.2]; Y1, Y2 = 7.4, 0.3

def card(x, y, n, title, body, fc='#f6f8fc'):
    """Rounded card: numbered title on top, icon area in the middle (drawn separately), body text at the bottom."""
    col = STEP_COL[n - 1]
    ax.scatter([x + 0.42], [y + H - 0.45], s=150, color=col, lw=0, zorder=3)
    ax.text(x + 0.42, y + H - 0.45, str(n), ha='center', va='center', color='white', fontsize=6.2, fontweight='bold', zorder=4)
    ax.text(x + 0.85, y + H - 0.45, title, va='center', fontsize=6.4, fontweight='bold', color=col)
    ax.text(x + 0.3, y + 0.25, body, va='bottom', fontsize=5.0, color=INK, linespacing=1.3)

def arr(p, q):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=10, color=LIGHT, lw=1.2, zorder=0))

# 1 isoform set
card(X[0], Y1, 1, 'Human isoform set',
     'UniProtKB/Swiss-Prot, reviewed\n- named gene, at least 20 residues\n- protein evidence (PE < 5)\n- genes with at least 2 isoforms\n\n'
     f'{CC["isoforms"]:,} isoforms, {CC["genes"]:,} genes')
cx, cy = X[0] + 2.3, Y1 + 4.0
for k in range(3):
    ax.add_patch(Rectangle((cx - 0.8, cy - 0.5 + k * 0.4), 1.6, 0.4, fc='#dfe8f6', ec=BLUE, lw=0.7, zorder=2))
    ax.add_patch(Ellipse((cx, cy - 0.1 + k * 0.4), 1.6, 0.3, fc='#dfe8f6', ec=BLUE, lw=0.7, zorder=2))
# 2 annotate
card(X[1], Y1, 2, 'Annotate genes',
     'HumanTFs; UniProt keywords,\nfamilies, EC numbers, location;\nMobiDB-lite, DisProt, CD-CODE\n\n'
     f'11 protein classes\n{PC["I_n_families_ge15"]} families with at least 15 genes')
for k, c in enumerate([BLUE, AQUA, ORANGE]):
    x0 = X[1] + 1.35 + k * 0.75; y0 = Y1 + 3.55
    ax.add_patch(Polygon([[x0, y0], [x0 + 0.55, y0], [x0 + 0.55, y0 + 0.8], [x0 + 0.35, y0 + 1.0], [x0, y0 + 1.0]], fc='white', ec=c, lw=0.9, zorder=2))
    for j in range(3): ax.plot([x0 + 0.1, x0 + 0.45], [y0 + 0.2 + j * 0.2] * 2, color=c, lw=0.5)
# 3 disorder
card(X[2], Y1, 3, 'Predict disorder',
     'MetaPredict v3.1, per residue\nscore of 0.5 or more: disordered\nIDR: at least 20 disordered\nresidues in a row\n\n'
     f'{CC["idrs"]:,} IDRs')
t = np.linspace(0, 1, 120); sc = 0.5 + 0.42 * np.sin(2 * np.pi * (t * 1.6 + 0.15)) * np.exp(-((t - 0.55) ** 2) / 0.5)
xx = X[2] + 0.5 + t * 3.2; yy = Y1 + 3.4 + sc * 1.3
ax.fill_between(xx, Y1 + 3.4, yy, where=sc >= 0.5, color=PINK, lw=0, zorder=2); ax.plot(xx, yy, color=PINK_E, lw=0.9, zorder=3)
ax.plot([xx[0], xx[-1]], [Y1 + 4.05] * 2, color=MUTED, lw=0.5, ls='--'); ax.text(xx[-1] + 0.08, Y1 + 4.05, '0.5', fontsize=4.6, va='center', color=MUTED)
# 4 transcripts with attrition bars
card(X[3], Y1, 4, 'Match transcripts',
     'candidate Ensembl transcripts per\nisoform; keep exact GENCODE v50\ntranslations that resolve to one\ncoding structure\n\nbars: isoforms kept at each stage')
steps = [('cohort', CC['isoforms']), ('with a transcript', CC['with_candidate_enst']), ('exact sequence', CC['exact_match']), ('one coding structure', CC['main'])]
for k, (lab, n) in enumerate(steps):
    yb = Y1 + 4.95 - k * 0.42; wb = 1.2 * n / CC['isoforms']
    ax.add_patch(Rectangle((X[3] + 0.3, yb), wb, 0.3, fc=BLUE, alpha=0.3 + 0.2 * k, lw=0, zorder=2))
    ax.text(X[3] + 0.38 + wb, yb + 0.15, f'{n:,} {lab}', fontsize=4.7, va='center', zorder=3)
for k in range(3): arr((X[k] + W + 0.05, Y1 + H / 2), (X[k + 1] - 0.05, Y1 + H / 2))
arr((X[3] + W / 2, Y1 - 0.05), (X[3] + W / 2, Y2 + H + 0.05))
# 5 project
card(X[3], Y2, 5, 'Project to GRCh38',
     'exon and coding coordinates;\nIDR and ordered segments\nprojected base by base, then\ntranslated back to protein\n\n'
     f'{CC["segments"]:,} of {CC["segments"]:,} segments pass')
seq = 'ATGGCTCAGCAT'
for k, ch in enumerate(seq):
    ax.add_patch(Rectangle((X[3] + 0.45 + k * 0.29, Y2 + 4.05), 0.27, 0.4, fc=GREEN if k < 6 else PINK, ec='white', lw=0.4, zorder=2))
    ax.text(X[3] + 0.585 + k * 0.29, Y2 + 4.25, ch, fontsize=4.4, ha='center', va='center', zorder=3)
for k, aa in enumerate('MAQH'):
    ax.add_patch(Rectangle((X[3] + 0.45 + k * 0.87, Y2 + 3.4), 0.83, 0.35, fc='white', ec=GREEN_E if k < 2 else PINK_E, lw=0.6, zorder=2))
    ax.text(X[3] + 0.865 + k * 0.87, Y2 + 3.575, aa, fontsize=5.0, ha='center', va='center')
ax.plot([X[3] + 4.0, X[3] + 4.12, X[3] + 4.35], [Y2 + 3.6, Y2 + 3.45, Y2 + 3.85], color='#2a9d8f', lw=1.6, solid_capstyle='round')
# 6 usage and positions
card(X[2], Y2, 6, 'Exons and positions',
     'C: exon in every mapped transcript\nA: in only some; residues matched\nby genomic codon: shared or\nisoform-specific\n\n'
     f'{CC["analysis_unique_exons"]:,} exons\n{RF["positions"] / 1e6:.1f} M coding positions')
ax.plot([X[2] + 0.4, X[2] + 4.2], [Y2 + 4.1] * 2, color=MUTED, lw=0.7, zorder=1)
for k, alt in enumerate([0, 1, 0, 1]):
    ax.add_patch(FancyBboxPatch((X[2] + 0.55 + k * 0.92, Y2 + 3.91), 0.72, 0.38, boxstyle='round,pad=0,rounding_size=0.06',
                                fc=PINK if alt else GREEN, ec=ORANGE if alt else INK, lw=0.9 if alt else 0.6, zorder=2))
    ax.text(X[2] + 0.91 + k * 0.92, Y2 + 4.1, 'A' if alt else 'C', fontsize=5.2, ha='center', va='center', color=ORANGE if alt else INK, zorder=3)
# 7 dating
card(X[1], Y2, 7, 'Gene age, splice sites',
     'gene age: GenOrigin\nsplice-site conservation: GT/AG\nkept in 405 aligned mammals\n\n'
     f'{EV["E1_genes_with_age"]:,} genes dated\n{EA["internal_with_both_sites"]:,} internal exons scored')
tx, ty = X[1] + 1.2, Y2 + 3.45
for (a, b), (c, d) in [((0, 0.5), (0.6, 0.5)), ((0.6, 0.15), (0.6, 0.85)), ((0.6, 0.85), (1.3, 0.85)), ((0.6, 0.15), (1.0, 0.15)), ((1.0, 0.0), (1.0, 0.3)),
                       ((1.0, 0.0), (1.6, 0.0)), ((1.0, 0.3), (1.6, 0.3)), ((1.3, 0.7), (1.3, 1.0)), ((1.3, 0.7), (1.6, 0.7)), ((1.3, 1.0), (1.6, 1.0))]:
    ax.plot([tx + a * 1.2, tx + c * 1.2], [ty + b, ty + d], color=MUTED, lw=0.8)
# 8 analysis
card(X[0], Y2, 8, 'Analyse by gene',
     'genes as the unit: resampling genes,\nfixed effects, clustered errors,\npaired tests; BH correction\n\nby class, family, gene age\nand splice-site conservation')
rng = np.random.default_rng(3); px = rng.uniform(0, 1, 18); py = 0.25 + 0.6 * px + rng.normal(0, 0.08, 18)
ax.plot([X[0] + 0.9, X[0] + 3.9], [Y2 + 3.35] * 2, color=MUTED, lw=0.6); ax.plot([X[0] + 0.9] * 2, [Y2 + 3.35, Y2 + 4.75], color=MUTED, lw=0.6)
ax.scatter(X[0] + 0.95 + px * 2.9, Y2 + 3.35 + py * 1.3, s=5, color=ORANGE, zorder=3)
ax.plot([X[0] + 1.0, X[0] + 3.8], [Y2 + 3.35 + 0.39, Y2 + 3.35 + 1.14], color=ORANGE, lw=0.7, ls='--')
for k in (3, 2, 1): arr((X[k] - 0.05, Y2 + H / 2), (X[k - 1] + W + 0.05, Y2 + H / 2))
save(fig, 'extended_data/ed1_methods')
