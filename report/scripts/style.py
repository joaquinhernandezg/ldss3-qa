"""Shared matplotlib style and paths for the report figures."""
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPORT = Path(__file__).resolve().parents[1]
FIGS = REPORT / 'figures'
TABLES = REPORT / 'tables'
RESULTS = REPORT / 'results'
REDUX = Path(os.environ.get('LDSS3_REDUX',
                            '/panoramix1/estudiantes/jhernandez/LDSS3_REDUCTION/report_2026/redux'))
RAW = Path(os.environ.get('LDSS3_RAW',
                          '/panoramix1/estudiantes/jhernandez/LDSS3_REDUCTION/report_2026/raw'))

# Categorical colours, always assigned in this order
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948')
INK, INK2, GRID = '#0b0b0b', '#52514e', '#d9d8d4'

AMP_COLOR = {1: BLUE, 2: ORANGE}
ION_COLOR = {'HeI': BLUE, 'NeI': ORANGE, 'ArI': AQUA, 'OH': VIOLET}

# Datasets analysed in the report: name -> (grism, mode, label)
DATASETS = {
    'vph_blue_longslit': ('VPH-Blue', 'longslit', 'VPH-Blue, 1.0" longslit'),
    'vph_red_longslit': ('VPH-Red', 'longslit', 'VPH-Red, 1.0" longslit'),
    'vph_all_longslit': ('VPH-All', 'longslit', 'VPH-All, 1.0" longslit'),
    'vph_blue_mos': ('VPH-Blue', 'mos', 'VPH-Blue, mask ATM3a2_1'),
    'vph_red_mos': ('VPH-Red', 'mos', 'VPH-Red, mask CDFS01'),
}

# Text width of the report (article, 1in margins on A4) in inches
TEXTWIDTH = 6.3


def setup():
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Computer Modern Roman', 'DejaVu Serif'],
        'mathtext.fontset': 'cm',
        'font.size': 8.5,
        'axes.labelsize': 8.5,
        'axes.titlesize': 8.5,
        'legend.fontsize': 7.5,
        'xtick.labelsize': 7.5,
        'ytick.labelsize': 7.5,
        'axes.edgecolor': INK2,
        'axes.labelcolor': INK,
        'xtick.color': INK2,
        'ytick.color': INK2,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'ytick.right': True,
        'axes.linewidth': 0.6,
        'lines.linewidth': 0.8,
        'axes.grid': False,
        'legend.frameon': False,
        'savefig.dpi': 200,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
        'figure.facecolor': 'white',
        'image.origin': 'lower',
    })
    FIGS.mkdir(exist_ok=True)
    TABLES.mkdir(exist_ok=True)


def save(fig, name):
    fig.savefig(FIGS / f'{name}.pdf')
    plt.close(fig)
    print(f'wrote {FIGS / name}.pdf')
