#!/usr/bin/env python3
"""Compact summary and h->0 extrapolation diagnostics for full_constants.py.

Reads agent_artifacts/full_constants_n*.json and writes
results/full_constants/full_constants_summary.json: per-mesh constants and checks
and least-squares fits of L+a/n^2, L+a/n^2+b/n^3 and L+a/n+b/n^2.
The fits are extrapolation diagnostics, not error bounds.
"""
import json
from pathlib import Path

import numpy as np

from paths import AGENT_ARTIFACTS

RESULTS = Path(__file__).parent/'results'/'full_constants'
MESHES = (5, 10, 15, 20, 30)
MODELS = {
    'L+a/n^2': lambda n: [np.ones_like(n), n**-2.],
    'L+a/n^2+b/n^3': lambda n: [np.ones_like(n), n**-2., n**-3.],
    'L+a/n+b/n^2': lambda n: [np.ones_like(n), n**-1., n**-2.],
}
PROVENANCE = dict(
    model='claude-opus-5-5',
    git_head='8c56dfdaeff1ff530c2d3598d77db868c6bd62db (branch johnson-mercier; '
             'full_constants.py and this script uncommitted)',
    software='Python 3.14.7, numpy 2.4.6, scipy 1.16.2 (SuperLU), FlexiBLAS',
    machine='AMD Ryzen 7 5825U, 22 GB RAM, Linux 7.1.13-200.fc44',
    reproduce=['cd addisons_experiments/bdm_hx_decomposition',
               'python3 full_constants.py --mesh N --refinements 8 16 32   # N=5,10,15',
               'python3 full_constants.py --mesh 20 --refinements 8 16',
               'systemd-run --user --scope -p MemoryMax=15G -p MemorySwapMax=1G '
               'python3 full_constants.py --mesh 30 --refinements 8 16 --chunk 128',
               'python3 fit_full_constants.py'],
    notes='lu_nnz (n<=20) is L.nnz+U.nnz; lu_storage_nnz (n=30) is SuperLU storage nnz.')
SUBSETS = {'all': MESHES, 'n>=10': (10, 15, 20, 30), 'n>=15': (15, 20, 30)}


def fit(points, model):
    n = np.array([p[0] for p in points], float)
    y = np.array([p[1] for p in points])
    basis = np.array(MODELS[model](n)).T
    if len(n) < basis.shape[1]:
        return None
    coef, *_ = np.linalg.lstsq(basis, y, rcond=None)
    return dict(L=float(coef[0]), coefficients=[float(c) for c in coef[1:]],
                max_abs_residual=float(np.abs(basis @ coef - y).max()),
                points=[[int(a), float(b)] for a, b in zip(n, y)])


def main():
    records = {n: json.loads((AGENT_ARTIFACTS/f'full_constants_n{n}.json').read_text())
               for n in MESHES if (AGENT_ARTIFACTS/f'full_constants_n{n}.json').exists()}
    meshes = {}
    for n, d in records.items():
        c = d['CoptHX']
        m = dict(bdm_dimension=d['bdm_dimension'], CoptHX=c['CoptHX'],
                 CoptHX_checks={k: c[k] for k in ('rayleigh', 'split_energy_ratio',
                                                  'reconstruction_relative',
                                                  'eigen_residual_relative', 'binv_symmetry')},
                 CoptHX_split_components=c['split_components'],
                 C_nolocal_discrete=c['C_nolocal_discrete'], refinements={})
        for r in d['refinements']:
            keep = ('spline_dimension', 'Cavg', 'Cavg_rayleigh', 'Cavg_eigen_residual_relative',
                    'Cavg_components_at_maximizer', 'Cavg_mean_subtracted',
                    'PiQ_mean_operator_norm', 'smooth_part_operator_norm_squared',
                    'potential_part_operator_norm_squared', 'patch_only_averaged',
                    'reconstruction_max_error', 'dominance_min_eig_Favg_minus_Binv',
                    'stationarity_residual', 'energy_identity_block_error', 'CoptHX_le_Cavg',
                    'Ccont', 'compared_to_refinement', 'product_map_change_norm',
                    'map_change_components', 'previous_Cavg_maximizer_ratio_now',
                    'continuous_map_change_norm', 'continuous_energies_at_Cavg_maximizer',
                    'lu_nnz', 'lu_storage_nnz', 'elapsed_seconds', 'peak_rss_gb')
            m['refinements'][str(r['refinement'])] = {k: r[k] for k in keep if k in r}
        meshes[str(n)] = m
    series = {'CoptHX': [(n, meshes[str(n)]['CoptHX']) for n in records]}
    for r in (8, 16):
        pts = [(n, meshes[str(n)]['refinements'][str(r)]['Cavg']) for n in records
               if str(r) in meshes[str(n)]['refinements']]
        series[f'Cavg_r{r}'] = pts
    best = []
    for n in records:
        rs = meshes[str(n)]['refinements']
        rmax = max(rs, key=int)
        best.append((n, rs[rmax]['Cavg'], int(rmax)))
    series['Cavg_best_r'] = [(n, y) for n, y, _ in best]
    fits = {}
    for name, pts in series.items():
        for sub, allowed in SUBSETS.items():
            chosen = [p for p in pts if p[0] in allowed]
            for model in MODELS:
                res = fit(chosen, model)
                if res is not None:
                    fits.setdefault(name, {}).setdefault(sub, {})[model] = res
    out = dict(description='Full squared product-norm constants, alpha=1, unit square, '
                           'uniform BDM1 triangulation with n x n squares. Fits are '
                           'extrapolation diagnostics only.',
               provenance=PROVENANCE,
               best_refinement_used=[[n, r] for n, _, r in best],
               meshes=meshes, series={k: [[int(a), float(b)] for a, b in v] for k, v in series.items()},
               fits=fits)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS/'full_constants_summary.json').write_text(json.dumps(out, indent=1)+'\n')
    for name, bysub in fits.items():
        for sub, bymodel in bysub.items():
            print(name, sub, {k: round(v['L'], 6) for k, v in bymodel.items()})


if __name__ == '__main__':
    main()
