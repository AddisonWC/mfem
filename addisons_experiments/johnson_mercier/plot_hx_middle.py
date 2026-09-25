#!/usr/bin/env python3
"""Plot selected smoother-centered and hybrid spectra from the retained CSV."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent
with (root / 'results/hx_middle.csv').open() as f:
    rows = list(csv.DictReader(f))
with (root / 'results/hx_spectrum.csv').open() as f:
    old = list(csv.DictReader(f))
fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
for ax, (patch, h1) in zip(axes, [('large', 'unsplit'), ('small', 'unsplit'), ('small', 'split')]):
    s = 1/3 if patch == 'large' else (.4 if h1 == 'unsplit' else .6)
    for order, c, color in [('HCSCH', 1, 'tab:blue'), ('CHSHC', 1, 'tab:cyan'), ('H(C+S)H', .5, 'tab:orange')]:
        selected = {}
        for row in rows:
            if (row['patch'] == patch and row['h1'] == h1 and row['composition'] == order
                and row['converged'] == '1' and abs(float(row['wS'])-s)<1e-9
                and float(row['wH']) == 1 and float(row['wC']) == c):
                r = int(row['refinement'])
                # Prefer tighter requested residual tolerance, then later validation.
                if r not in selected or float(row['tolerance']) <= float(selected[r]['tolerance']):
                    selected[r] = row
        pairs = sorted((r, 1/float(row['lambda_min']) if patch=='large' and order!='H(C+S)H' else float(row['kappa'])) for r, row in selected.items())
        ax.plot([4*2**r for r,k in pairs], [k for r,k in pairs], 'o-', color=color, label=order)
    ref = sorted((int(row['refinement']),float(row['kappa'])) for row in old
                 if row['patch']==patch and row['h1']==h1 and row['config']=='chosen_mult'
                 and row['converged']=='1' and row['seed']=='17')
    ax.plot([4*2**r for r,k in ref], [k for r,k in ref], 's--', color='tab:green', label='SHCHS reference')
    ax.set_title(f'{patch} patches, {h1} H1')
    ax.set_xscale('log',base=2)
    ax.set_xlabel('Square subdivisions per side')
    ax.set_ylabel('Spectral condition number')
    ax.grid(alpha=.25)
axes[0].legend(fontsize=8)
fig.savefig(root / 'results/hx_middle.png', dpi=180)
