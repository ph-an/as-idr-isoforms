"""Shared figure style for the manuscript figures."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]; FIG = REPO / 'manuscript' / 'figures'; AN = REPO / 'analysis'
BLUE, ORANGE, AQUA, MUTED, GRID, INK, LIGHT, RED = '#2a78d6', '#eb6834', '#1baf7a', '#52514e', '#e6e5e0', '#0b0b0b', '#c9c8c2', '#c2410c'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.6, 'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.labelcolor': INK})
SHORT = {'Metabolic enzymes': 'Metabolic enzymes', 'Transporters': 'Transporters', 'Secreted': 'Secreted', 'Kinases & phosphatases': 'Kinases & phosphatases',
         'Receptors & channels': 'Receptors & channels', 'Other membrane': 'Other membrane', 'Cytoskeletal & structural': 'Cytoskeletal & structural',
         'Other / unclassified': 'Other / unclassified', 'Chromatin & RNA processing': 'Chromatin & RNA proc.',
         'Transcription co-regulators': 'Transcription co-regulators', 'DNA-binding TFs': 'DNA-binding TFs'}
AGE = ['<100', '100-500', '500-1000', '>1000']

def tag(ax, t, x=-0.12, y=1.04):
    ax.text(x, y, t, transform=ax.transAxes, fontweight='bold', fontsize=9, va='bottom')

def grid(ax, axis='y'):
    ax.grid(axis=axis, color=GRID, lw=0.5); ax.set_axisbelow(True)

def save(fig, name):
    """Write manuscript/figures/<name>.pdf and .png; name includes the subfolder, e.g. 'main/fig1_overview'."""
    (FIG / name).parent.mkdir(parents=True, exist_ok=True)
    for e in ['pdf', 'png']:
        fig.savefig(FIG / f'{name}.{e}', dpi=300, bbox_inches='tight')
    print('written', FIG / f'{name}.pdf')

REF = '#a8a7a1'   # reference lines (zero, expectation, equality): light grey, dashed, always labelled

def refline(ax, v, axis='x', label=None, ha=None, fs=5, y=1.01, x=0.99, va='bottom'):
    """Dashed light-grey reference line with a short label saying what it marks.
    axis='x': vertical line at x=v, label above the axes; axis='y': horizontal line at y=v, label at the right end."""
    if axis == 'x':
        ax.axvline(v, color=REF, lw=0.7, ls='--', zorder=0.5)
        if label: ax.text(v, y, label, transform=ax.get_xaxis_transform(), ha=ha or 'center', va='bottom', fontsize=fs, color=MUTED)
    else:
        ax.axhline(v, color=REF, lw=0.7, ls='--', zorder=0.5)
        if label: ax.text(x, v, label, transform=ax.get_yaxis_transform(), ha=ha or 'right', va=va, fontsize=fs, color=MUTED)

# colours used by the schematics in Fig. 1a-b and Extended Data Fig. 1 (ordered vs disordered sequence, as on the poster)
ORD, ORD_E, DIS, DIS_E = '#bfe0a8', '#5f9e43', '#f4a6a6', '#d9534f'
