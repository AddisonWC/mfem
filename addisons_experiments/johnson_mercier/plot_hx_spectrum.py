"""Plot the selected square-mesh HX spectral study (requires pandas/matplotlib)."""
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data', type=Path, default=HERE/'results/hx_spectrum.csv')
parser.add_argument('--output', type=Path, default=HERE/'output/hx_spectrum.png')
args = parser.parse_args()
df = pd.read_csv(args.data)
df = df[(df.converged == 1) & (df.seed == 17) & (df.kappa > 0)]
styles = {
    'base': ('Additive (1, 1, 1)', '#52525b', 'o'),
    's033': ('Additive (.33, 1, 1)', '#d97706', 's'),
    'chosen_add': ('Additive: large (.75,2,1); small unsplit (1,1,1), split (.75,1,1)', '#7c3aed', '^'),
    'mult_SHCHS': ('SHCHS (.33 large / .5 small, 1, 1)', '#0284c7', 'D'),
    'chosen_mult': ('SHCHS (.4 large / .6 small, 1, 1)', '#059669', 'v'),
}
fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True, sharey="row")
for row, patch in enumerate(['large', 'small']):
    for col, h1 in enumerate(['unsplit', 'split']):
        ax = axes[row, col]
        sub = df[(df.patch == patch) & (df.h1 == h1)]
        for cfg, (label, color, marker) in styles.items():
            g = sub[sub.config == cfg].drop_duplicates('refinement', keep='last').sort_values('refinement')
            ax.plot(4*2**g.refinement, g.kappa, color=color, marker=marker,
                    label=label, linewidth=1.6, markersize=4)
        ax.set_title(f'{patch.capitalize()} patches, {h1} H1')
        ax.set_xscale('log', base=2)
        ax.set_xticks([4,8,16,32,64,128], labels=['4','8','16','32','64','128'])
        ax.set_ylabel(r'$\kappa(BA)$')
        ax.grid(alpha=.22)
        ax.set_ylim(bottom=0)
        if row == 1:
            ax.set_xlabel('Square subdivisions per side (h = 1/n)')
handles, labels = axes[0,0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False, fontsize=9)
fig.suptitle('Johnson–Mercier HX: square-mesh refinement', fontsize=15)
fig.tight_layout(rect=(0, .10, 1, .96))
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=180)
print(args.output)
