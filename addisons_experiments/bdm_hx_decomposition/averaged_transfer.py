"""Area-weighted averaging of trianglewise P1 L2 projections.

For a vector field s, (J s)(v) = sum_{K containing v} int_K (12 lam_v - 3) s
divided by |omega_v|.  On each K, (12 lam_v - 3)/|K| is the L2-dual of lam_v in
P1, so J reproduces continuous P1 exactly.  Applied componentwise.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse

import plot_decompositions as base
import sparse_continuous


def subtriangle_rule(refinement, points=5):
    """Collapsed Gauss rule on the r^2 uniform subtriangles of the reference."""
    z, wz = leggauss(points)
    z, wz = (z+1)/2, wz/2
    a, b = np.meshgrid(z, z, indexing='ij')
    wa, wb = np.meshgrid(wz, wz, indexing='ij')
    ref = np.stack(((1-a)*(1-b), a, (1-a)*b), axis=-1).reshape(-1, 3)
    weights = (wa*wb*(1-a)).ravel()
    r = refinement
    bary = lambda i, j: np.array([1-(i+j)/r, i/r, j/r])
    subs = []
    for i in range(r):
        for j in range(r-i):
            a0, b0, c0 = bary(i, j), bary(i+1, j), bary(i, j+1)
            subs.append(np.stack((a0, b0, c0)))
            if i+j < r-1:
                subs.append(np.stack((b0, bary(i+1, j+1), c0)))
    lam = np.concatenate([ref @ sub for sub in subs])
    return lam, np.tile(weights, r*r)/(r*r)


def patch_areas(xy, tri):
    _, area = base.geometry(xy, tri)
    return np.bincount(tri.ravel(), np.repeat(area, 3), minlength=len(xy))


def averaged_rotgrad(potential, xy, tri, refinement):
    """Sparse map from spline coefficients to interleaved J(rotgrad psi).

    Subtriangles of size h/r lie in single spline cells, and the integrand
    has degree at most 6 there, so the degree-9 collapsed rule is exact.
    """
    _, area = base.geometry(xy, tri)
    lam, weights = subtriangle_rule(refinement)
    dual = 12*lam-3
    rows, blocks = [], []
    for k, t in enumerate(tri):
        _, gx, gy, *_ = sparse_continuous.evaluate(potential, lam @ xy[t])
        w = (2*area[k]*weights)[:, None]*dual
        blocks += [sparse.csr_matrix(w.T) @ gy, -(sparse.csr_matrix(w.T) @ gx)]
        rows += [2*t, 2*t+1]
    order = np.concatenate(rows)
    stacked = sparse.vstack(blocks).tocsr()
    select = sparse.csr_matrix((np.ones(len(order)), (order, np.arange(len(order)))),
                               shape=(2*len(xy), len(order)))
    scale = sparse.diags(np.repeat(1/patch_areas(xy, tri), 2))
    return (scale @ select @ stacked).tocsr()


def averaged_p1(xy, tri):
    """J applied to interleaved continuous vector P1 coefficients (exact rule)."""
    gradlam, area = base.geometry(xy, tri)
    lam, weights = subtriangle_rule(1)
    nv = len(xy)
    j = np.zeros((2*nv, 2*nv))
    for k, t in enumerate(tri):
        local = (lam*(2*area[k]*weights)[:, None]).T @ (12*lam-3)
        for c in range(2):
            j[np.ix_(2*t+c, 2*t+c)] += local.T
    return j/np.repeat(patch_areas(xy, tri), 2)[:, None]
