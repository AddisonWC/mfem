#!/usr/bin/env python3
"""Nodal versus area-averaged smooth transfer after the same continuous solve."""
import argparse
import json
import time

import numpy as np
from scipy import linalg

import averaged_transfer
import maximize_patch_ratio as search
from paths import AGENT_ARTIFACTS
from refine_continuous import largest


def top(a, b):
    val, vec = linalg.eigh((a+a.T)/2, b, subset_by_index=(len(b)-1, len(b)-1),
                           check_finite=False)
    u = vec[:, 0]
    return float(val[0]), (u if u[np.argmax(np.abs(u))] > 0 else -u)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mesh', type=int, default=5)
    parser.add_argument('--refinements', nargs='+', type=int, default=[4, 8, 16])
    args = parser.parse_args()
    report = AGENT_ARTIFACTS/f'averaging_comparison_n{args.mesh}.json'
    records = json.loads(report.read_text()) if report.exists() else []
    previous = None
    for r in args.refinements:
        start = time.monotonic()
        res = search.find_worst(args.mesh, basis_type='spline', spline_refinement=r)
        xy, tri, pre, emap = res['xy'], res['tri'], res['pre'], res['emap']
        a, w, sd, qd = res['hdiv'].toarray(), res['wlocal'], pre['sd'], pre['qd']
        potential = search.TensorSplinePotential(args.mesh*r)
        jmap = averaged_transfer.averaged_rotgrad(potential, xy, tri, r)
        jp1 = averaged_transfer.averaged_p1(xy, tri)
        maps = dict(nodal=res['localmap'],
                    averaged=qd @ (res['scalar'] @ emap) - sd @ (jmap @ emap))
        # Smooth coefficient maps u -> s0 and u -> interpolated s.
        lagrange = linalg.solve(pre['bs']+pre['bq'], np.eye(len(a)), assume_a='sym')
        s0 = linalg.cho_solve(pre['sf'], sd.T @ lagrange)
        q0 = linalg.cho_solve(pre['qf'], qd.T @ lagrange)
        smooth = dict(nodal=s0+res['nodal'] @ emap, averaged=s0+jmap @ emap)
        hs = pre['hs']
        record = dict(mesh=args.mesh, refinement=r, spline_dimension=len(emap),
                      stationarity_residual=res['residual'],
                      p1_reproduction_error=float(np.abs(jp1-np.eye(len(jp1))).max()))
        for name, lmap in maps.items():
            ratio, u = top(lmap.T @ w @ lmap, a)
            local = lmap @ u
            psi = emap @ u
            s0u, q0u = s0 @ u, q0 @ u
            jpsi = (jmap if name == 'averaged' else res['nodal']) @ psi
            recon = sd @ (s0u+jpsi) + qd @ (q0u-res['scalar'] @ psi) + local
            # Uncancelled form Pi_V s - T s with T the transfer applied to all of s.
            ts0 = jp1 @ s0u if name == 'averaged' else s0u
            full = sd @ s0u + qd @ (res['scalar'] @ psi) - sd @ (ts0+jpsi)
            other = maps['nodal' if name == 'averaged' else 'averaged'] @ u
            info = dict(maximum_ratio=ratio,
                        rayleigh_check=float(local @ w @ local/(u @ a @ u)),
                        reconstruction_relative=float(np.linalg.norm(recon-u)/np.linalg.norm(u)),
                        uncancelled_local_relative=float(np.linalg.norm(full-local)/max(1e-300, np.linalg.norm(local))),
                        maximizer_ratio_under_other_transfer=float(other @ w @ other/(u @ a @ u)),
                        interpolated_s_h1_squared=float(smooth[name] @ u @ hs @ (smooth[name] @ u)),
                        smooth_transfer_h1_operator_norm_squared=largest(smooth[name].T @ hs @ smooth[name], a))
            if previous is not None:
                delta = lmap-previous['maps'][name]
                info['patch_operator_change_norm'] = float(np.sqrt(max(0, largest(delta.T @ w @ delta, a))))
                pu = previous['fields'][name]
                pl = lmap @ pu
                info['previous_maximizer_ratio_now'] = float(pl @ w @ pl/(pu @ a @ pu))
            assert abs(info['rayleigh_check']-ratio) < 1e-8*max(1, ratio)
            assert info['reconstruction_relative'] < 1e-8
            record[name] = info
            np.save(AGENT_ARTIFACTS/f'averaging_{name}_maximizer_n{args.mesh}_ref{r}_bdm_dofs.npy', u)
            record.setdefault('_fields', {})[name] = u
        diff = maps['averaged']-maps['nodal']
        record['nodal_maximizer_continuous_norms'] = res['continuous_norms']
        record['nodal_minus_averaged_operator_norm'] = float(np.sqrt(max(0, largest(diff.T @ w @ diff, a))))
        if previous is not None:
            record['compared_to_refinement'] = previous['refinement']
        record['elapsed_seconds'] = time.monotonic()-start
        fields = record.pop('_fields')
        records = [x for x in records if x['refinement'] != r]+[record]
        records.sort(key=lambda x: x['refinement'])
        report.write_text(json.dumps(records, indent=2)+'\n')
        print(json.dumps(record), flush=True)
        previous = dict(maps=maps, fields=fields, refinement=r)
        del res, maps, emap


if __name__ == '__main__':
    main()
