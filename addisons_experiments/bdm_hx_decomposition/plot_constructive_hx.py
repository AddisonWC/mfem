#!/usr/bin/env python3
"""Constructive HX splits using a commuting moment interpolation pair."""

import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import legder, legval, leggauss
from scipy import linalg

import plot_decompositions as base


def eval_p2_gradient(fine, points):
    xy, tri, edges, te, phi = fine
    n = round(np.sqrt(len(xy))) - 1
    pts = np.asarray(points).reshape(-1, 2)
    x = np.clip(pts[:, 0], 0, 1-1e-13)
    y = np.clip(pts[:, 1], 0, 1-1e-13)
    i = np.minimum((x*n).astype(int), n-1)
    j = np.minimum((y*n).astype(int), n-1)
    fx, fy = x*n-i, y*n-j
    tid = 2*(j*n+i)+(fy>fx).astype(int)
    gradlam, _ = base.geometry(xy, tri)
    out = np.empty((len(pts), 2))
    for k in np.unique(tid):
        mask = tid == k
        p = xy[tri[k]]
        lam12 = (pts[mask]-p[0]) @ np.linalg.inv(
            np.stack((p[1]-p[0], p[2]-p[0]), axis=-1)).T
        lam = np.column_stack((1-lam12.sum(axis=1), lam12))
        g = gradlam[k]
        vals = phi[np.r_[tri[k], len(xy)+te[k]]]
        deriv = np.zeros((len(lam), 2))
        for a in range(3):
            deriv += vals[a]*(4*lam[:, a]-1)[:, None]*g[a]
        for m, (a, b) in enumerate(((0, 1), (1, 2), (2, 0)), 3):
            deriv += 4*vals[m]*(lam[:, a, None]*g[b] + lam[:, b, None]*g[a])
        out[mask] = deriv
    return out.reshape(np.shape(points))


def eval_bdm(u, points, n, td, basis):
    pts = np.asarray(points).reshape(-1, 2)
    x = np.clip(pts[:, 0], 0, 1-1e-13)
    y = np.clip(pts[:, 1], 0, 1-1e-13)
    i = np.minimum((x*n).astype(int), n-1)
    j = np.minimum((y*n).astype(int), n-1)
    tid = 2*(j*n+i) + (y*n-j > x*n-i).astype(int)
    c = np.einsum("tij,tj->ti", basis, u[td])[tid]
    return np.column_stack((c[:, 0]+x*c[:, 1]+y*c[:, 2],
                            c[:, 3]+x*c[:, 4]+y*c[:, 5])).reshape(np.shape(points))


class SmoothPotential:
    def __init__(self, degree=6):
        self.modes = [(i, j) for total in range(1, degree+1)
                      for i in range(total+1) for j in [total-i]]
        self.coef = np.zeros(len(self.modes))

    def basis(self, points):
        pts = np.asarray(points).reshape(-1, 2)
        x, y = 2*pts[:, 0]-1, 2*pts[:, 1]-1
        order = max(i+j for i, j in self.modes)
        polys_x, polys_y, deriv_x, deriv_y = [], [], [], []
        for k in range(order+1):
            c = np.eye(order+1)[k]
            polys_x.append(legval(x, c)); polys_y.append(legval(y, c))
            deriv_x.append(2*legval(x, legder(c)))
            deriv_y.append(2*legval(y, legder(c)))
        val = np.column_stack([polys_x[i]*polys_y[j] for i, j in self.modes])
        gradx = np.column_stack([deriv_x[i]*polys_y[j] for i, j in self.modes])
        grady = np.column_stack([polys_x[i]*deriv_y[j] for i, j in self.modes])
        return val, np.stack((gradx, grady), axis=1)

    def value_grad(self, points):
        val, grad = self.basis(points)
        return val @ self.coef, np.einsum("ncm,m->nc", grad, self.coef)

    def basis_hessian(self, points):
        pts = np.asarray(points).reshape(-1, 2)
        x, y = 2*pts[:, 0]-1, 2*pts[:, 1]-1
        order = max(i+j for i, j in self.modes)
        px, py, dx, dy, dxx, dyy = [], [], [], [], [], []
        for k in range(order+1):
            c = np.eye(order+1)[k]
            px.append(legval(x, c)); py.append(legval(y, c))
            dx.append(2*legval(x, legder(c)))
            dy.append(2*legval(y, legder(c)))
            dxx.append(4*legval(x, legder(c, 2)))
            dyy.append(4*legval(y, legder(c, 2)))
        val = np.column_stack([px[i]*py[j] for i, j in self.modes])
        gx = np.column_stack([dx[i]*py[j] for i, j in self.modes])
        gy = np.column_stack([px[i]*dy[j] for i, j in self.modes])
        hxx = np.column_stack([dxx[i]*py[j] for i, j in self.modes])
        hxy = np.column_stack([dx[i]*dy[j] for i, j in self.modes])
        hyy = np.column_stack([px[i]*dyy[j] for i, j in self.modes])
        return val, gx, gy, hxx, hxy, hyy


def no_local_coefficients(u, pre):
    b = pre["bs"] + pre["bq"]
    lagrange = linalg.solve(b, u, assume_a="sym", check_finite=False)
    s = linalg.cho_solve(pre["sf"], pre["sd"].T @ lagrange,
                          check_finite=False)
    q = linalg.cho_solve(pre["qf"], pre["qd"].T @ lagrange,
                          check_finite=False)
    assert np.linalg.norm(u - pre["sd"]@s - pre["qd"]@q) < 1e-9
    return s, q


def fit_to_helmholtz(u, discrete_s, fine, tri, xy, td, basis):
    # L2 fit of rotgrad ψ to (Neumann gradient - discrete no-local smooth).
    points = np.concatenate([base.QPTS @ xy[t] for t in tri])
    target = eval_p2_gradient(fine, points) - eval_bdm(discrete_s, points, 10, td, basis)
    psi = SmoothPotential(8)
    _, grads = psi.basis(points)
    rot = np.stack((grads[:, 1], -grads[:, 0]), axis=1)
    matrix = rot.reshape(-1, len(psi.modes))
    # A small degree penalty suppresses high-frequency overfitting to the
    # fine-grid Galerkin gradient's elementwise jumps.
    penalty = np.array([(i+j)**2 for i, j in psi.modes], float)
    matrix = np.vstack((matrix, 2e-2*np.diag(penalty)))
    rhs = np.r_[target.ravel(), np.zeros(len(psi.modes))]
    psi.coef = np.linalg.lstsq(matrix, rhs, rcond=None)[0]
    fit = rot.reshape(-1, len(psi.modes)) @ psi.coef
    err = float(np.linalg.norm(target.ravel()-fit) / max(np.linalg.norm(target), 1e-14))
    return psi, err


def minimize_continuous_h1(scoef, qcoef, xy, tri, te, degree=8):
    """Galerkin minimizer of ||s0+R∇ψ||²_H1 + ||q0-ψ||²_H1."""
    psi = SmoothPotential(degree)
    m = len(psi.modes)
    gram = np.zeros((m, m))
    linear = np.zeros(m)
    constant = 0.0
    gram_s = np.zeros((m, m))
    linear_s = np.zeros(m)
    constant_s = 0.0
    nv = len(xy)
    gradlam, area = base.geometry(xy, tri)
    z, wz = leggauss(degree+2)
    z, wz = (z+1)/2, wz/2
    aa, bb = np.meshgrid(z, z, indexing="ij")
    wa, wb = np.meshgrid(wz, wz, indexing="ij")
    a, b = aa.ravel(), bb.ravel()
    ref_lam = np.column_stack((1-a-(1-a)*b, a, (1-a)*b))
    weights_ref = (wa*wb*(1-aa)).ravel()
    for k, t in enumerate(tri):
        points = ref_lam @ xy[t]
        weights = 2*area[k]*weights_ref
        val, gx, gy, hxx, hxy, hyy = psi.basis_hessian(points)
        # Nine terms: two vector values, four vector derivatives, scalar
        # value, and two scalar derivatives.
        modes = np.stack((gy, -gx, hxy, hyy, -hxx, -hxy,
                          -val, -gx, -gy), axis=1)
        svertex = scoef.reshape(-1, 2)[t]
        svalues = ref_lam @ svertex
        ds = svertex.T @ gradlam[k]
        qids = np.r_[t, nv+te[k]]
        qlocal = qcoef[qids]
        qvalues = np.empty(len(points))
        dq = np.empty((len(points), 2))
        for j, lam in enumerate(ref_lam):
            shape, dshape = base.p2shape(lam, gradlam[k])
            qvalues[j] = shape @ qlocal
            dq[j] = qlocal @ dshape
        target = np.column_stack((svalues,
            np.tile(ds.ravel(), (len(points), 1)), qvalues, dq))
        gram += np.einsum("qim,qin,q->mn", modes, modes, weights)
        linear += np.einsum("qim,qi,q->m", modes, target, weights)
        constant += np.einsum("qi,qi,q->", target, target, weights)
        gram_s += np.einsum("qim,qin,q->mn", modes[:, :6], modes[:, :6], weights)
        linear_s += np.einsum("qim,qi,q->m", modes[:, :6], target[:, :6], weights)
        constant_s += np.einsum("qi,qi,q->", target[:, :6], target[:, :6], weights)
    psi.coef = linalg.solve(gram, -linear, assume_a="pos", check_finite=False)
    minimum = float(constant + 2*linear@psi.coef + psi.coef@gram@psi.coef)
    smooth_energy = float(constant_s + 2*linear_s@psi.coef +
                          psi.coef@gram_s@psi.coef)
    return psi, minimum, smooth_energy, minimum-smooth_energy


def p2_moment_interpolant(psi, xy, edges):
    # The P2 edge DOF is its edge average. Store the equivalent midpoint value.
    nv = len(xy)
    vertex_val, _ = psi.value_grad(xy)
    t, w = leggauss(7)
    t, w = (t+1)/2, w/2
    midvals = np.empty(len(edges))
    for e, (a, b) in enumerate(edges):
        points = xy[a][None, :] + t[:, None]*(xy[b]-xy[a])[None, :]
        vals, _ = psi.value_grad(points)
        avg = w @ vals
        midvals[e] = 1.5*avg - .25*(vertex_val[a]+vertex_val[b])
    return np.r_[vertex_val, midvals]


def bdm_moment_interpolant(scoef, psi, xy, edges, normal):
    # Moment projection on each edge: integrate normal flux against P1.
    t, w = leggauss(7)
    t, w = (t+1)/2, w/2
    nodal = scoef.reshape(-1, 2)
    out = np.empty(2*len(edges))
    gram_inv = np.array([[4., -2.], [-2., 4.]])
    for e, (a, b) in enumerate(edges):
        points = xy[a][None, :] + t[:, None]*(xy[b]-xy[a])[None, :]
        _, grad = psi.value_grad(points)
        rot = np.column_stack((grad[:, 1], -grad[:, 0]))
        smooth = (1-t[:, None])*nodal[a] + t[:, None]*nodal[b] + rot
        flux = smooth @ normal[e]
        moments = np.array([np.dot(w, flux*(1-t)), np.dot(w, flux*t)])
        out[2*e:2*e+2] = gram_inv @ moments
    return out


def construct(u, pre, xy, tri, te, edges, normal, td, basis, hdiv,
              continuous_mode="helmholtz_fit"):
    scoef, qcoef = no_local_coefficients(u, pre)
    discrete_s = pre["sd"] @ scoef
    if continuous_mode == "helmholtz_fit":
        fine = base.neumann_helmholtz(u, 10, td, basis, 40)
        psi, fit_error = fit_to_helmholtz(u, discrete_s, fine, tri, xy, td, basis)
        continuous_energy = None
        convergence = None
        continuous_s_energy = continuous_q_energy = None
    elif continuous_mode == "h1_optimal":
        convergence = {}
        for degree in (4, 6, 8):
            psi, continuous_energy, continuous_s_energy, continuous_q_energy = minimize_continuous_h1(
                scoef, qcoef, xy, tri, te, degree)
            convergence[str(degree)] = continuous_energy
        fit_error = None
    else:
        raise ValueError(continuous_mode)
    psi_q = p2_moment_interpolant(psi, xy, edges)
    pi_v_s = bdm_moment_interpolant(scoef, psi, xy, edges, normal)
    commuting_defect = float(np.linalg.norm(pi_v_s - discrete_s - pre["qd"]@psi_q))
    assert commuting_defect < 1e-9, commuting_defect
    _, vertex_grad = psi.value_grad(xy)
    nodal_s = scoef.reshape(-1, 2) + np.column_stack((vertex_grad[:, 1],
                                                       -vertex_grad[:, 0]))
    nodal_s = nodal_s.ravel()
    smooth = pre["sd"] @ nodal_s
    q_h = qcoef - psi_q
    solenoidal = pre["qd"] @ q_h
    tilde = pi_v_s - smooth
    reconstruction = float(np.linalg.norm(u - smooth - solenoidal - tilde))
    assert reconstruction < 1e-9, reconstruction
    patches = []
    local_norm2 = 0.0
    for ids, _ in pre["fact"]:
        val = tilde[ids]
        patches.append((ids, val))
        local_norm2 += float(val @ pre["av"][np.ix_(ids, ids)] @ val)
    smooth_norm2 = float(nodal_s @ pre["hs"] @ nodal_s)
    q_norm2 = float(q_h @ pre["hq"] @ q_h)
    original = float(u @ hdiv @ u)
    metrics = dict(
        source_hdiv_squared=original,
        smooth_h1_squared=smooth_norm2,
        potential_h1_squared=q_norm2,
        local_hdiv_squared=local_norm2,
        total_alpha_quarter=smooth_norm2+q_norm2+.25*local_norm2,
        total_alpha_one=smooth_norm2+q_norm2+local_norm2,
        total_alpha_four=smooth_norm2+q_norm2+4*local_norm2,
        reconstruction_l2_dof=reconstruction,
        commuting_l2_dof=commuting_defect)
    if continuous_mode == "helmholtz_fit":
        metrics["relative_helmholtz_fit_l2"] = fit_error
    else:
        metrics["continuous_h1_product_squared"] = continuous_energy
        metrics["continuous_s_h1_squared"] = continuous_s_energy
        metrics["continuous_q_h1_squared"] = continuous_q_energy
        metrics["nodal_s_h1_squared_ratio"] = smooth_norm2/continuous_s_energy
        metrics["moment_q_h1_squared_ratio"] = q_norm2/continuous_q_energy
        metrics["moment_s_hdiv_squared_ratio"] = (
            float(pi_v_s @ hdiv @ pi_v_s)/continuous_s_energy)
        metrics["polynomial_degree_convergence"] = convergence
    return smooth, solenoidal, patches, metrics


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    xy, tri, edges, te, _ = base.mesh(10)
    normal, td, basis = base.bdm_data(xy, tri, edges, te)
    hdiv = base.bdm_hdiv(xy, tri, td, basis)
    h1 = base.assemble_h1(xy, tri, te, len(xy), 1)
    h2 = base.assemble_h1(xy, tri, te, len(xy), 2)
    smap, qmap = base.prolongations(xy, edges, normal)
    pre = base.preconditioner(smap, qmap, h1, h2, hdiv, edges, len(xy))
    result = {}
    result_h1 = {}
    for name, u in base.source_fields(xy, edges, normal, qmap).items():
        _, _, _, _, _, optimum = base.optimal_split(u, pre, 1.)
        source_im, tid, xx, yy = base.evaluate_bdm_image(u, 10, td, basis)
        for mode in ("helmholtz_fit", "h1_optimal"):
            smooth, sol, patches, metrics = construct(
                u, pre, xy, tri, te, edges, normal, td, basis, hdiv, mode)
            metrics["optimal_alpha_one_squared"] = optimum["weighted_total"]
            metrics["constructive_over_optimal_alpha_one"] = (
                metrics["total_alpha_one"] / optimum["weighted_total"])
            smooth_im, _, _, _ = base.evaluate_bdm_image(smooth, 10, td, basis)
            sol_im, _, _, _ = base.evaluate_bdm_image(sol, 10, td, basis)
            local_im = base.patch_energy_image(patches, td, basis, tid, xx, yy)
            info = ("Moment BDM + P2 interpolants\n"
                    f"Source H(div) norm²: {metrics['source_hdiv_squared']:.4g}\n"
                    f"P1 H¹ norm²: {metrics['smooth_h1_squared']:.4g}\n"
                    f"P2 H¹ norm²: {metrics['potential_h1_squared']:.4g}\n"
                    f"Σ patch H(div) norm²: {metrics['local_hdiv_squared']:.4g}\n"
                    f"Total, α=¼ / 1 / 4: {metrics['total_alpha_quarter']:.3g} / "
                    f"{metrics['total_alpha_one']:.3g} / {metrics['total_alpha_four']:.3g}\n"
                    f"α=1 squared ratio to optimum: "
                    f"{metrics['constructive_over_optimal_alpha_one']:.2f}×\n"
                    f"Reconstruction: {metrics['reconstruction_l2_dof']:.1e}\n")
            if mode == "helmholtz_fit":
                info += f"Helmholtz fit error: {metrics['relative_helmholtz_fit_l2']:.1%}"
                path = out / f"{name}_constructive_hx.png"
                title = f"{name.replace('_', ' ').title()} · constructive HX split"
                result[name] = metrics
            else:
                c = metrics["polynomial_degree_convergence"]
                info += (f"Continuous product norm²: {c['8']:.4g}\n"
                         f"Degrees 4/6/8: {c['4']:.3g}/{c['6']:.3g}/{c['8']:.3g}\n"
                         f"ΠS / s H¹ norm²: "
                         f"{metrics['nodal_s_h1_squared_ratio']:.2f}×")
                path = out / f"{name}_constructive_h1.png"
                title = f"{name.replace('_', ' ').title()} · H¹-optimal continuous split"
                result_h1[name] = metrics
            base.plot_figure(path, title, source_im, smooth_im, sol_im,
                             local_im, info)
            print(name, mode, json.dumps(metrics), flush=True)
    (out / "constructive_metrics.json").write_text(json.dumps(result, indent=2)+"\n")
    (out / "constructive_h1_metrics.json").write_text(
        json.dumps(result_h1, indent=2)+"\n")


if __name__ == "__main__":
    main()
