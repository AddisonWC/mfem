#!/usr/bin/env python3
"""Plot completed continuous-refinement diagnostics; no extrapolated values."""
import json
from pathlib import Path
import matplotlib
from paths import AGENT_ARTIFACTS
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    out = AGENT_ARTIFACTS
    fig, axes = plt.subplots(1,2,figsize=(16,9),layout='constrained')
    for n in (5,10,15):
        path = out/f'continuous_refinement_n{n}.json'
        if not path.exists():
            continue
        data = json.loads(path.read_text())
        axes[0].plot([d['refinement'] for d in data],
                     [d['maximum_ratio'] for d in data], 'o-',label=f'{n}×{n} BDM mesh')
        changes = [d for d in data if 'patch_operator_change_norm' in d]
        axes[1].plot([d['refinement'] for d in changes],
                     [d['patch_operator_change_norm'] for d in changes], 'o-',label=f'{n}×{n} BDM mesh')
    for ax in axes:
        ax.set_xscale('log',base=2)
        ax.set_xticks([1,2,4,8,16,32],labels=['1','2','4','8','16','32'])
        ax.set_yscale('log')
        ax.set_xlabel('Spline cells per BDM grid cell')
        ax.grid(True,which='both',alpha=.25)
        ax.legend()
    axes[0].set_ylabel('Maximum squared patch norm / squared source norm')
    axes[0].set_title('Worst patch amplification')
    axes[1].set_ylabel('H(div) → patch product norm of map difference')
    axes[1].set_title('Change from preceding tested refinement')
    fig.suptitle('Continuous H¹-product solve: refinement still matters',fontsize=20)
    fig.savefig(out/'continuous_resolution.png',dpi=120)
    plt.close(fig)


if __name__ == '__main__':
    main()
