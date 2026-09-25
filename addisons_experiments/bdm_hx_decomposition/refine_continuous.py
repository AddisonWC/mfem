#!/usr/bin/env python3
"""Nested continuous resolution study, including operator changes and fixed inputs."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy import linalg
import maximize_patch_ratio as search
from paths import AGENT_ARTIFACTS


def largest(a, b):
    return float(linalg.eigh((a+a.T)/2, b, subset_by_index=(len(b)-1,len(b)-1),
                            eigvals_only=True, check_finite=False)[0])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mesh', type=int, default=10)
    parser.add_argument('--refinements', nargs='+', type=int, default=[1,2,4,8])
    args = parser.parse_args()
    out = AGENT_ARTIFACTS
    report = out/f'continuous_refinement_n{args.mesh}.json'
    records = json.loads(report.read_text()) if report.exists() else []
    previous = None
    for r in args.refinements:
        start = time.monotonic()
        result = search.find_worst(args.mesh, basis_type='spline', spline_refinement=r)
        a, w = result['hdiv'].toarray(), result['wlocal']
        pre = result['pre']
        b = pre['bs']+pre['bq']
        energy = linalg.inv(b)+result['rhs'].T @ result['emap']
        record = dict(mesh=args.mesh, refinement=r, spline_dimension=len(result['emap']),
                      maximum_ratio=result['ratio'], stationarity_residual=result['residual'])
        u = result['field']
        record['maximizer_continuous_energy'] = float(u @ energy @ u)
        if previous is not None:
            record['compared_to_refinement'] = previous['refinement']
            delta = result['localmap']-previous['localmap']
            record['patch_operator_change_norm'] = np.sqrt(max(0,largest(delta.T @ w @ delta,a)))
            record['continuous_energy_drop_operator_norm'] = largest(previous['energy']-energy,a)
            pu = previous['field']
            local = result['localmap'] @ pu
            record['previous_maximizer_ratio_now'] = float(local @ w @ local/(pu @ a @ pu))
        record['plot_metrics'] = search.plot_worst(result, out)
        norms = result['continuous_norms']
        assert np.isclose(sum(norms.values()), record['maximizer_continuous_energy'],
                          rtol=1e-7, atol=1e-8)
        assert record['plot_metrics']['reconstruction_l2_dof'] < 1e-7
        tag = f'n{args.mesh}_spline_ref{r}'
        np.save(out/f'worst_patch_ratio_{tag}_bdm_dofs.npy',u)
        record['elapsed_seconds'] = time.monotonic()-start
        records = [item for item in records if item['refinement'] != r]+[record]
        records.sort(key=lambda item: item['refinement'])
        report.write_text(json.dumps(records,indent=2)+'\n')
        print(json.dumps({k:v for k,v in record.items() if k!='plot_metrics'}),flush=True)
        previous = {k:result[k] for k in ('localmap','field')}
        previous['energy'] = energy
        previous['refinement'] = r
        del result


if __name__ == '__main__':
    main()
