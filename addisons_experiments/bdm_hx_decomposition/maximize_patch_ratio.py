#!/usr/bin/env python3
"""Largest vertex-patch amplification of the constructive H1 HX procedure."""

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import linalg
from scipy.interpolate import BSpline
from scipy.sparse.linalg import spsolve

import plot_decompositions as base
import plot_constructive_hx as constructive


class TensorSplinePotential:
    """C2 bicubic tensor splines on the square grid, a conforming H2 space."""

    def __init__(self, n):
        self.n = n
        self.knots = np.r_[np.zeros(4), np.arange(1, n)/n, np.ones(4)]
        self.nb = len(self.knots)-4
        self.spline = BSpline(self.knots, np.eye(self.nb), 3, axis=0)
        self.modes = range(self.nb*self.nb)

    def basis_hessian(self, points):
        pts = np.asarray(points).reshape(-1, 2)
        bx = [self.spline(pts[:, 0], nu=k) for k in range(3)]
        by = [self.spline(pts[:, 1], nu=k) for k in range(3)]
        pair = lambda a, b: np.einsum("ni,nj->nij", a, b).reshape(len(pts), -1)
        return (pair(bx[0], by[0]), pair(bx[1], by[0]),
                pair(bx[0], by[1]), pair(bx[2], by[0]),
                pair(bx[1], by[1]), pair(bx[0], by[2]))

    def basis(self, points):
        val, gx, gy, _, _, _ = self.basis_hessian(points)
        return val, np.stack((gx, gy), axis=1)


def continuous_quadratic_maps(xy, tri, te, degree, potential=None,
                              integration_refinement=1):
    """Build quadratic H1 minimization in a global polynomial H2 correction."""
    if potential is None:
        potential = constructive.SmoothPotential(degree)
    m, nv, ne = len(potential.modes), len(xy), int(te.max())+1
    gram = np.zeros((m, m))
    cs = np.zeros((m, 2*nv))
    cq = np.zeros((m, nv+ne))
    gradlam, area = base.geometry(xy, tri)
    z, wz = leggauss(degree+2)
    z, wz = (z+1)/2, wz/2
    aa, bb = np.meshgrid(z, z, indexing="ij")
    wa, wb = np.meshgrid(wz, wz, indexing="ij")
    a, b = aa.ravel(), bb.ravel()
    lam_ref = np.column_stack((1-a-(1-a)*b, a, (1-a)*b))
    wref = (wa*wb*(1-aa)).ravel()
    nq = len(lam_ref)
    r = integration_refinement
    bary = lambda i, j: np.array([1-(i+j)/r, i/r, j/r])
    subtriangles = []
    for i in range(r):
        for j in range(r-i):
            a0, b0, c0 = bary(i, j), bary(i+1, j), bary(i, j+1)
            subtriangles.append(np.stack((a0, b0, c0)))
            if i+j < r-1:
                subtriangles.append(np.stack((b0, bary(i+1, j+1), c0)))
    sub_lams = [lam_ref @ sub for sub in subtriangles]
    for k, t in enumerate(tri):
        sids = np.ravel(np.column_stack((2*t, 2*t+1)))
        qids = np.r_[t, nv+te[k]]
        for lam in sub_lams:
            x = lam @ xy[t]
            weight = 2*area[k]/(r*r)*wref
            val, gx, gy, hxx, hxy, hyy = potential.basis_hessian(x)
            modes = np.stack((gy, -gx, hxy, hyy, -hxx, -hxy,
                              -val, -gx, -gy), axis=1)
            active = np.flatnonzero(np.any(np.abs(modes) > 1e-13, axis=(0, 1)))
            modes = modes[:, :, active]
            gram[np.ix_(active, active)] += np.einsum(
                "qim,qin,q->mn", modes, modes, weight)
            sbasis = np.zeros((nq, 6, 6))
            for v in range(3):
                sbasis[:, 0, 2*v] = lam[:, v]
                sbasis[:, 1, 2*v+1] = lam[:, v]
                sbasis[:, 2, 2*v] = gradlam[k, v, 0]
                sbasis[:, 3, 2*v] = gradlam[k, v, 1]
                sbasis[:, 4, 2*v+1] = gradlam[k, v, 0]
                sbasis[:, 5, 2*v+1] = gradlam[k, v, 1]
            qbasis = np.zeros((nq, 3, 6))
            for j, barypoint in enumerate(lam):
                shape, grad = base.p2shape(barypoint, gradlam[k])
                qbasis[j, 0] = shape
                qbasis[j, 1:] = grad.T
            cs[np.ix_(active, sids)] += np.einsum(
                "qim,qij,q->mj", modes[:, :6], sbasis, weight)
            cq[np.ix_(active, qids)] += np.einsum(
                "qim,qij,q->mj", modes[:, 6:], qbasis, weight)
    return potential, gram, cs, cq


def moment_maps(potential, xy, edges, integration_refinement=1):
    """Maps polynomial coefficients to P2 moments and nodal P1 vector data."""
    nv, ne = len(xy), len(edges)
    m = len(potential.modes)
    vertices, grad = potential.basis(xy)
    nodal = np.empty((2*nv, m))
    nodal[0::2] = grad[:, 1]
    nodal[1::2] = -grad[:, 0]
    scalar = np.empty((nv+ne, m))
    scalar[:nv] = vertices
    t, w = leggauss(7)
    t, w = (t+1)/2, w/2
    r = integration_refinement
    t = ((np.arange(r)[:, None]+t[None, :])/r).ravel()
    w = np.tile(w/r, r)
    points = xy[edges[:, 0], None, :] + t[None, :, None] * (
        xy[edges[:, 1]]-xy[edges[:, 0]])[:, None, :]
    edgeval, _ = potential.basis(points.reshape(-1, 2))
    average = np.einsum("etm,t->em", edgeval.reshape(ne, len(t), m), w)
    scalar[nv:] = 1.5*average - .25*(vertices[edges[:, 0]]+
                                     vertices[edges[:, 1]])
    return nodal, scalar


def find_worst(n, degree=8, basis_type="polynomial", spline_refinement=1):
    xy, tri, edges, te, _ = base.mesh(n)
    normal, td, basis = base.bdm_data(xy, tri, edges, te)
    hdiv = base.bdm_hdiv(xy, tri, td, basis)
    h1 = base.assemble_h1(xy, tri, te, len(xy), 1)
    h2 = base.assemble_h1(xy, tri, te, len(xy), 2)
    smap, qmap = base.prolongations(xy, edges, normal)
    pre = base.preconditioner(smap, qmap, h1, h2, hdiv, edges, len(xy))
    correction_space = (TensorSplinePotential(n*spline_refinement)
                        if basis_type == "spline"
                        else None)
    potential, gram, cs, cq = continuous_quadratic_maps(
        xy, tri, te, degree, correction_space,
        spline_refinement if basis_type == "spline" else 1)
    nodal, scalar = moment_maps(potential, xy, edges,
                                 spline_refinement if basis_type == "spline" else 1)
    dmap = pre["qd"] @ scalar - pre["sd"] @ nodal
    # psi coefficients = E u, where u is in endpoint-normal BDM coordinates.
    b = pre["bs"] + pre["bq"]
    bfactor = linalg.cho_factor(b, check_finite=False)
    rhs = (cs @ linalg.cho_solve(pre["sf"], pre["sd"].T, check_finite=False) +
           cq @ linalg.cho_solve(pre["qf"], pre["qd"].T, check_finite=False))
    rhs = linalg.cho_solve(bfactor, rhs.T, check_finite=False).T
    emap = -linalg.cho_solve(linalg.cho_factor(gram, check_finite=False),
                              rhs, check_finite=False)
    # W is the sum of exact local H(div) matrices after endpoint partitioning.
    wlocal = np.zeros_like(pre["av"])
    for ids, _ in pre["fact"]:
        wlocal[np.ix_(ids, ids)] = pre["av"][np.ix_(ids, ids)]
    ksmall = dmap.T @ wlocal @ dmap
    ksmall = (ksmall+ksmall.T)/2
    eig, vec = linalg.eigh(ksmall, check_finite=False)
    active = eig > max(eig[-1]*1e-11, 1e-12)
    factor = vec[:, active] * np.sqrt(eig[active])[None, :]
    ainverse_et = spsolve(hdiv, emap.T)
    small = factor.T @ emap @ ainverse_et @ factor
    small = (small+small.T)/2
    top, z = linalg.eigh(small, subset_by_index=(small.shape[0]-1,
                                                 small.shape[0]-1),
                         check_finite=False)
    ratio = float(top[0])
    u = ainverse_et @ (factor @ z[:, 0]) / np.sqrt(ratio)
    if u[np.argmax(np.abs(u))] < 0:
        u = -u
    local = dmap @ (emap @ u)
    exact_ratio = float(local @ wlocal @ local / (u @ hdiv @ u))
    assert abs(exact_ratio-ratio) < 1e-8*max(1, ratio)
    return dict(n=n, degree=degree, basis_type=basis_type,
                spline_refinement=spline_refinement,
                ratio=ratio, field=u,
                xy=xy, tri=tri, edges=edges, te=te, normal=normal,
                td=td, basis=basis, hdiv=hdiv, pre=pre,
                dmap=dmap, emap=emap, nodal=nodal, scalar=scalar)


def plot_worst(result, out):
    n = result["n"]
    u = result["field"]
    if result["basis_type"] == "polynomial":
        smooth, sol, patches, metrics = constructive.construct(
            u, result["pre"], result["xy"], result["tri"], result["te"],
            result["edges"], result["normal"], result["td"], result["basis"],
            result["hdiv"], "h1_optimal", result["degree"])
    else:
        pre = result["pre"]
        scoef, qcoef = constructive.no_local_coefficients(u, pre)
        psi = result["emap"] @ u
        smooth_coef = scoef + result["nodal"] @ psi
        q_coef = qcoef - result["scalar"] @ psi
        smooth = pre["sd"] @ smooth_coef
        sol = pre["qd"] @ q_coef
        local = result["dmap"] @ psi
        patches = [(ids, local[ids]) for ids, _ in pre["fact"]]
        local_norm2 = sum(float(vals @ pre["av"][np.ix_(ids, ids)] @ vals)
                          for ids, vals in patches)
        metrics = dict(source_hdiv_squared=float(u @ result["hdiv"] @ u),
                       local_hdiv_squared=local_norm2,
                       smooth_h1_squared=float(smooth_coef @ pre["hs"] @ smooth_coef),
                       potential_h1_squared=float(q_coef @ pre["hq"] @ q_coef),
                       reconstruction_l2_dof=float(np.linalg.norm(u-smooth-sol-local)))
    numerator = metrics["local_hdiv_squared"]
    denominator = metrics["source_hdiv_squared"]
    assert abs(numerator/denominator-result["ratio"]) < 1e-7
    src, tid, xx, yy = base.evaluate_bdm_image(u, n, result["td"], result["basis"])
    sim, _, _, _ = base.evaluate_bdm_image(smooth, n, result["td"], result["basis"])
    qim, _, _, _ = base.evaluate_bdm_image(sol, n, result["td"], result["basis"])
    density = base.patch_energy_image(patches, result["td"], result["basis"],
                                      tid, xx, yy)
    info = (f"Worst generalized eigenvector\n{n}×{n} cells · BDM1\n"
            f"Σ patch H(div) norm²: {numerator:.4g}\n"
            f"Source H(div) norm²: {denominator:.4g}\n"
            f"Ratio: {result['ratio']:.4g}\n"
            f"P1 H¹ norm²: {metrics['smooth_h1_squared']:.4g}\n"
            f"P2 H¹ norm²: {metrics['potential_h1_squared']:.4g}\n"
            f"Reconstruction: {metrics['reconstruction_l2_dof']:.1e}\n"
            f"Continuous correction: {result['basis_type']} "
            f"{result['degree'] if result['basis_type']=='polynomial' else ''}")
    tag = (f"degree{result['degree']}" if result["basis_type"] == "polynomial"
           else f"spline_ref{result['spline_refinement']}")
    base.plot_figure(out / f"worst_patch_ratio_n{n}_{tag}.png",
                     f"Worst patch amplification · {n}×{n} mesh · {tag}",
                     src, sim, qim, density, info)
    return metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--meshes", nargs="+", type=int, default=[5, 10, 15])
    parser.add_argument("--plot-n", type=int, default=10)
    parser.add_argument("--degree", type=int, default=8)
    parser.add_argument("--basis", choices=("polynomial", "spline"),
                        default="polynomial")
    parser.add_argument("--spline-refinement", type=int, default=1)
    args = parser.parse_args()
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    summary = {}
    for n in args.meshes:
        result = find_worst(n, args.degree, args.basis,
                            args.spline_refinement)
        summary[str(n)] = {"mesh": n, "bdm_dimension": len(result["field"]),
                           "continuous_polynomial_degree": args.degree,
                           "correction_basis": args.basis,
                           "spline_refinement": args.spline_refinement,
                           "maximum_ratio": result["ratio"]}
        print(n, result["ratio"], flush=True)
        if n == args.plot_n:
            summary[str(n)]["plot_metrics"] = plot_worst(result, out)
            tag = (f"degree{args.degree}" if args.basis == "polynomial"
                   else f"spline_ref{args.spline_refinement}")
            np.save(out / f"worst_patch_ratio_n{n}_{tag}_bdm_dofs.npy",
                    result["field"])
    suffix = (args.basis if args.basis == "polynomial" else
              f"spline_ref{args.spline_refinement}")
    (out / f"worst_patch_ratio_{suffix}.json").write_text(
        json.dumps(summary, indent=2)+"\n")


if __name__ == "__main__":
    main()
