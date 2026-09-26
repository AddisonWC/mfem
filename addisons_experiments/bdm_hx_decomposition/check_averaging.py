#!/usr/bin/env python3
"""Exactness checks for the averaged P1 transfer on spline and P1 inputs."""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.interpolate import BSpline

import averaged_transfer as avg
import maximize_patch_ratio as search
import sparse_continuous
from paths import AGENT_ARTIFACTS


def spline_coefficients(potential, fx_list):
    """Tensor Greville collocation; exact for polynomials of degree <= 3."""
    k = potential.knots
    greville = (k[1:-3]+k[2:-2]+k[3:-1])/3
    col = BSpline.design_matrix(greville, k, 3).toarray()
    return sum(np.outer(np.linalg.solve(col, fx(greville)),
                        np.linalg.solve(col, fy(greville))).ravel()
               for fx, fy in fx_list)


def direct_average(xy, tri, grad):
    """Independent J(rotgrad psi) with one high-order rule per triangle."""
    gradlam, area = search.base.geometry(xy, tri)
    lam, weights = avg.subtriangle_rule(1, points=8)
    out = np.zeros(2*len(xy))
    for k, t in enumerate(tri):
        x = lam @ xy[t]
        gx, gy = grad(x[:, 0], x[:, 1])
        wd = (2*area[k]*weights)[:, None]*(12*lam-3)
        out[2*t] += wd.T @ gy
        out[2*t+1] -= wd.T @ gx
    return out/np.repeat(avg.patch_areas(xy, tri), 2)


def main():
    results = {}
    for n, r in ((3, 2), (5, 4)):
        xy, tri, edges, te, _ = search.base.mesh(n)
        p = search.TensorSplinePotential(n*r)
        jmap = avg.averaged_rotgrad(p, xy, tri, r)
        nodal, _ = sparse_continuous.moments(p, xy, edges, r)
        one = lambda x: np.ones_like(x)
        # psi = 1 + x + 2y + .3x^2 - .7xy + .5y^2: rotgrad psi is affine.
        quad = [(one, one), (lambda x: x, one), (one, lambda y: 2*y),
                (lambda x: .3*x**2, one), (lambda x: -.7*x, lambda y: y),
                (one, lambda y: .5*y**2)]
        c = spline_coefficients(p, quad)
        affine = float(np.abs(jmap @ c - nodal @ c).max())
        # psi = x^3 - 2x^2 y + x y^2 + .5 y^3: rotgrad psi is quadratic.
        cubic = [(lambda x: x**3, one), (lambda x: -2*x**2, lambda y: y),
                 (lambda x: x, lambda y: y**2), (one, lambda y: .5*y**3)]
        c = spline_coefficients(p, cubic)
        grad = lambda x, y: (3*x**2-4*x*y+y**2, -2*x**2+2*x*y+1.5*y**2)
        ref = direct_average(xy, tri, grad)
        quadratic = float(np.abs(jmap @ c - ref).max()/np.abs(ref).max())
        # Generic splines: doubling the knot-aligned subdivision changes nothing.
        coeff = np.random.default_rng(7).standard_normal(p.nb**2)
        fine = avg.averaged_rotgrad(p, xy, tri, 2*r) @ coeff
        subdivision = float(np.abs(jmap @ coeff - fine).max()/np.abs(fine).max())
        p1 = float(np.abs(avg.averaged_p1(xy, tri)-np.eye(2*len(xy))).max())
        item = dict(p1_reproduction=p1, affine_rotgrad_vs_nodal=affine,
                    quadratic_rotgrad_vs_direct_relative=quadratic,
                    random_spline_subdivision_relative=subdivision)
        for key, value in item.items():
            assert value < 1e-11, (n, r, key, value)
        results[f'n{n}_ref{r}'] = item
        print(n, r, item)
    (AGENT_ARTIFACTS/'averaging_checks.json').write_text(json.dumps(results, indent=2)+'\n')


if __name__ == '__main__':
    main()
