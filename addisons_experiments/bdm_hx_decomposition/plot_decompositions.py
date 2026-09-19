#!/usr/bin/env python3
"""Optimal BDM1 HX splits and an approximate continuous Neumann Helmholtz split.

The BDM1 coordinates are normal values at edge endpoints, with one globally
fixed normal per edge. This basis makes the vertex patch spaces a direct sum.
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import spsolve


QPTS = np.array([[1 / 6, 1 / 6, 2 / 3], [1 / 6, 2 / 3, 1 / 6],
                 [2 / 3, 1 / 6, 1 / 6]])
QWTS = np.full(3, 1 / 3)


def mesh(n):
    xy = np.array([(i / n, j / n) for j in range(n + 1)
                   for i in range(n + 1)], dtype=float)
    vertex = lambda i, j: j * (n + 1) + i
    tri = []
    for j in range(n):
        for i in range(n):
            a, b, c, d = vertex(i, j), vertex(i + 1, j), \
                         vertex(i + 1, j + 1), vertex(i, j + 1)
            tri.extend(((a, b, c), (a, c, d)))
    tri = np.asarray(tri, dtype=int)
    edge_ids = {}
    edges = []
    te = []
    for t in tri:
        local = []
        for x, y in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            pair = tuple(sorted((x, y)))
            if pair not in edge_ids:
                edge_ids[pair] = len(edges)
                edges.append(pair)
            local.append(edge_ids[pair])
        te.append(local)
    edges = np.asarray(edges, dtype=int)
    te = np.asarray(te, dtype=int)
    mid = (xy[edges[:, 0]] + xy[edges[:, 1]]) / 2
    return xy, tri, edges, te, mid


def geometry(xy, tri):
    p = xy[tri]
    mat = np.stack((p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=-1)
    inv = np.linalg.inv(mat)
    g1 = inv[:, 0, :]
    g2 = inv[:, 1, :]
    gradlam = np.stack((-g1 - g2, g1, g2), axis=1)
    area = np.linalg.det(mat) / 2
    return gradlam, area


def p2shape(lam, gradlam):
    vals = np.empty(6)
    grads = np.empty((6, 2))
    for a in range(3):
        vals[a] = lam[a] * (2 * lam[a] - 1)
        grads[a] = (4 * lam[a] - 1) * gradlam[a]
    for k, (a, b) in enumerate(((0, 1), (1, 2), (2, 0)), start=3):
        vals[k] = 4 * lam[a] * lam[b]
        grads[k] = 4 * (lam[a] * gradlam[b] + lam[b] * gradlam[a])
    return vals, grads


def assemble_h1(xy, tri, te, n_vertices, degree):
    gradlam, area = geometry(xy, tri)
    nd = n_vertices if degree == 1 else n_vertices + int(te.max()) + 1
    rows, cols, data = [], [], []
    for k, t in enumerate(tri):
        ids = t if degree == 1 else np.r_[t, n_vertices + te[k]]
        block = np.zeros((len(ids), len(ids)))
        for w, lam in zip(QWTS, QPTS):
            if degree == 1:
                vals, grads = lam, gradlam[k]
            else:
                vals, grads = p2shape(lam, gradlam[k])
            block += area[k] * w * (np.outer(vals, vals) + grads @ grads.T)
        ii, jj = np.meshgrid(ids, ids, indexing="ij")
        rows.extend(ii.ravel()); cols.extend(jj.ravel()); data.extend(block.ravel())
    return sparse.coo_matrix((data, (rows, cols)), shape=(nd, nd)).tocsr()


def bdm_data(xy, tri, edges, te):
    ne = len(edges)
    ev = xy[edges[:, 1]] - xy[edges[:, 0]]
    normal = np.column_stack((ev[:, 1], -ev[:, 0])) / np.linalg.norm(ev, axis=1)[:, None]
    tri_dofs = np.empty((len(tri), 6), dtype=int)
    basis = np.empty((len(tri), 6, 6))
    for k, t in enumerate(tri):
        vand = np.empty((6, 6))
        ids = []
        for e in te[k]:
            for side in range(2):
                v = edges[e, side]
                x, y = xy[v]
                nx, ny = normal[e]
                vand[len(ids)] = [nx, nx*x, nx*y, ny, ny*x, ny*y]
                ids.append(2*e + side)
        tri_dofs[k] = ids
        basis[k] = np.linalg.inv(vand)
    return normal, tri_dofs, basis


def bdm_hdiv(xy, tri, tri_dofs, basis):
    _, area = geometry(xy, tri)
    nd = int(tri_dofs.max()) + 1
    rows, cols, data = [], [], []
    for k, t in enumerate(tri):
        block = np.zeros((6, 6))
        div = basis[k, 1] + basis[k, 5]
        for w, lam in zip(QWTS, QPTS):
            x, y = lam @ xy[t]
            v = np.stack((basis[k, 0] + x*basis[k, 1] + y*basis[k, 2],
                          basis[k, 3] + x*basis[k, 4] + y*basis[k, 5]))
            block += area[k] * w * (v.T @ v + np.outer(div, div))
        ids = tri_dofs[k]
        ii, jj = np.meshgrid(ids, ids, indexing="ij")
        rows.extend(ii.ravel()); cols.extend(jj.ravel()); data.extend(block.ravel())
    return sparse.coo_matrix((data, (rows, cols)), shape=(nd, nd)).tocsr()


def prolongations(xy, edges, normal):
    nv, ne = len(xy), len(edges)
    rs, cs, ds, rq, cq, dq = [], [], [], [], [], []
    for e, (a, b) in enumerate(edges):
        length = np.linalg.norm(xy[b] - xy[a])
        for side, v in enumerate((a, b)):
            row = 2*e + side
            for c in range(2):
                rs.append(row); cs.append(2*v+c); ds.append(normal[e, c])
            for col, val in zip((a, nv+e, b),
                                ((-3, 4, -1) if side == 0 else (1, -4, 3))):
                rq.append(row); cq.append(col); dq.append(val / length)
    s = sparse.coo_matrix((ds, (rs, cs)), shape=(2*ne, 2*nv)).tocsr()
    q = sparse.coo_matrix((dq, (rq, cq)), shape=(2*ne, nv+ne)).tocsr()
    return s, q


def interpolate_field(xy, edges, normal, fn):
    out = np.empty(2*len(edges))
    for e, (a, b) in enumerate(edges):
        out[2*e] = normal[e] @ fn(*xy[a])
        out[2*e+1] = normal[e] @ fn(*xy[b])
    return out


def source_fields(xy, edges, normal, qmap, return_components=False):
    pi = np.pi
    def wave(x, y):
        return np.array([np.sin(pi*x)*np.cos(pi*y) + .25*x,
                         np.cos(pi*x)*np.sin(pi*y) - .15*y])
    def shear(x, y):
        return np.array([-.55*np.sin(2*pi*y)*(.7+.3*x) + .2*y,
                         .55*np.sin(2*pi*x)*(.7+.3*y) + .3*x])
    u1 = interpolate_field(xy, edges, normal, wave)
    u2 = interpolate_field(xy, edges, normal, shear)
    mid = (xy[edges[:, 0]] + xy[edges[:, 1]]) / 2
    # A localized P2 edge-midpoint pattern produces structured tangential jumps.
    envelope = np.exp(-25*((mid[:, 0]-.5)**2 + (mid[:, 1]-.5)**2))
    q = np.r_[np.zeros(len(xy)), .018 * envelope *
              np.cos(5*pi*mid[:, 0]) * np.sin(4*pi*mid[:, 1])]
    u3 = .35*u1 + qmap @ q
    # Four seeded Fourier bands are interpolated into continuous P2, then
    # normalized in the assembled H1 norm to match a smooth P1 vector field.
    _, tri, check_edges, te, _ = mesh(10)
    assert np.array_equal(edges, check_edges)
    h1 = assemble_h1(xy, tri, te, len(xy), 1)
    h2 = assemble_h1(xy, tri, te, len(xy), 2)
    smooth_values = np.array([.7*wave(*p) + .3*shear(*p) for p in xy])
    smooth_coef = smooth_values.ravel()
    smooth_h1 = sparse.kron(h1, sparse.eye(2, format="csr"), format="csr")
    smooth_norm = np.sqrt(smooth_coef @ smooth_h1 @ smooth_coef)
    nodes = np.vstack((xy, mid))
    rng = np.random.default_rng(20260919)
    q_noise = np.zeros(len(nodes))
    for lo, hi, weight in ((1, 2, .65), (3, 4, .48),
                           (5, 6, .34), (7, 9, .24)):
        band = np.zeros(len(nodes))
        for _ in range(8):
            kx, ky = (rng.integers(lo, hi+1, size=2) *
                      rng.choice((-1, 1), size=2))
            phase = rng.uniform(0, 2*pi)
            sign = rng.choice((-1., 1.))
            band += sign*np.cos(2*pi*(kx*nodes[:, 0] + ky*nodes[:, 1])+phase) \
                    / np.hypot(kx, ky)
        band -= band.mean()
        band /= np.sqrt(band @ h2 @ band)
        q_noise += weight*band
    q_noise -= q_noise.mean()
    q_noise *= smooth_norm / np.sqrt(q_noise @ h2 @ q_noise)
    u4 = qmap @ q_noise + prolongations(xy, edges, normal)[0] @ smooth_coef
    fields = {"compression_wave": u1, "shear_wave": u2,
              "localized_edge_pattern": u3,
              "multiscale_noisy_potential": u4}
    if return_components:
        return fields, dict(noisy_smooth_coefficients=smooth_coef,
                            noisy_potential_coefficients=q_noise,
                            smooth_h1_norm=smooth_norm,
                            potential_h1_norm=float(np.sqrt(q_noise @ h2 @ q_noise)))
    return fields


def preconditioner(smap, qmap, h1scalar, h1p2, hdiv, edges, nv):
    hs = sparse.kron(h1scalar, sparse.eye(2, format="csr"), format="csr").toarray()
    hq = h1p2.toarray()
    sd, qd = smap.toarray(), qmap.toarray()
    sf = linalg.cho_factor(hs, check_finite=False)
    qf = linalg.cho_factor(hq, check_finite=False)
    bs = sd @ linalg.cho_solve(sf, sd.T, check_finite=False)
    bq = qd @ linalg.cho_solve(qf, qd.T, check_finite=False)
    av = hdiv.toarray()
    bl = np.zeros_like(av)
    vertex_dofs = [[] for _ in range(nv)]
    for e, pair in enumerate(edges):
        for side, v in enumerate(pair):
            vertex_dofs[v].append(2*e+side)
    fact = []
    for ids in vertex_dofs:
        ids = np.asarray(ids, dtype=int)
        fac = linalg.cho_factor(av[np.ix_(ids, ids)], check_finite=False)
        bl[np.ix_(ids, ids)] = linalg.cho_solve(fac, np.eye(len(ids)), check_finite=False)
        fact.append((ids, fac))
    return dict(hs=hs, hq=hq, sd=sd, qd=qd, sf=sf, qf=qf,
                bs=bs, bq=bq, bl=bl, fact=fact, av=av)


def optimal_split(u, pre, alpha):
    b = pre["bs"] + pre["bq"]
    if alpha is not None:
        b = b + pre["bl"] / alpha
    lagrange = linalg.solve(b, u, assume_a="sym", check_finite=False)
    scoef = linalg.cho_solve(pre["sf"], pre["sd"].T @ lagrange,
                              check_finite=False)
    qcoef = linalg.cho_solve(pre["qf"], pre["qd"].T @ lagrange,
                              check_finite=False)
    smooth = pre["sd"] @ scoef
    solenoidal = pre["qd"] @ qcoef
    local = np.zeros_like(u)
    patches = None
    if alpha is not None:
        patches = []
        for ids, fac in pre["fact"]:
            coef = linalg.cho_solve(fac, lagrange[ids], check_finite=False) / alpha
            local[ids] = coef
            patches.append((ids, coef))
    residual = np.linalg.norm(u - smooth - solenoidal - local)
    smooth_energy = float(scoef @ pre["hs"] @ scoef)
    potential_energy = float(qcoef @ pre["hq"] @ qcoef)
    local_energy = 0.0
    if patches is not None:
        for ids, vals in patches:
            local_energy += float(vals @ pre["av"][np.ix_(ids, ids)] @ vals)
    product_energy = float(u @ lagrange)
    assert abs(product_energy - (smooth_energy + potential_energy +
                (0 if alpha is None else alpha*local_energy))) < 1e-8
    energies = dict(smooth_h1_squared=smooth_energy,
                    potential_h1_squared=potential_energy,
                    local_hdiv_squared=local_energy,
                    weighted_total=product_energy)
    return smooth, solenoidal, local, patches, residual, energies


def evaluate_bdm_image(coef, n, tri_dofs, basis, resolution=600):
    coords = (np.arange(resolution) + .5) / resolution
    xx, yy = np.meshgrid(coords, coords)
    ii = np.minimum((xx*n).astype(int), n-1)
    jj = np.minimum((yy*n).astype(int), n-1)
    fx, fy = xx*n-ii, yy*n-jj
    tid = 2*(jj*n+ii) + (fy > fx).astype(int)
    coeff = np.einsum("tij,tj->ti", basis, coef[tri_dofs])
    cc = coeff[tid]
    vx = cc[..., 0] + xx*cc[..., 1] + yy*cc[..., 2]
    vy = cc[..., 3] + xx*cc[..., 4] + yy*cc[..., 5]
    return np.stack((vx, vy)), tid, xx, yy


def patch_energy_image(patches, tri_dofs, basis, tid, xx, yy):
    energy = np.zeros_like(xx)
    if patches is None:
        return energy
    for ids, vals in patches:
        global_coef = np.zeros(int(tri_dofs.max())+1)
        global_coef[ids] = vals
        local = global_coef[tri_dofs]
        active = np.any(local != 0, axis=1)
        if not np.any(active):
            continue
        coeff = np.einsum("tij,tj->ti", basis[active], local[active])
        mask = active[tid]
        # Only triangles adjacent to this vertex are visited.
        indices = np.searchsorted(np.flatnonzero(active), tid[mask])
        cc = coeff[indices]
        x, y = xx[mask], yy[mask]
        vx = cc[:, 0] + x*cc[:, 1] + y*cc[:, 2]
        vy = cc[:, 3] + x*cc[:, 4] + y*cc[:, 5]
        div = cc[:, 1] + cc[:, 5]
        energy[mask] += vx*vx + vy*vy + div*div
    return energy


def neumann_helmholtz(u, n, coarse_dofs, coarse_basis, fine_n=40):
    """P2 Galerkin approximation to the continuous Neumann H1 projection."""
    xy, tri, edges, te, _ = mesh(fine_n)
    nv = len(xy)
    gradlam, area = geometry(xy, tri)
    rows, cols, data = [], [], []
    rhs = np.zeros(nv+len(edges))
    for k, t in enumerate(tri):
        ids = np.r_[t, nv+te[k]]
        block = np.zeros((6, 6))
        load = np.zeros(6)
        for w, lam in zip(QWTS, QPTS):
            x, y = lam @ xy[t]
            _, grads = p2shape(lam, gradlam[k])
            ci, cj = min(int(x*n), n-1), min(int(y*n), n-1)
            tid = 2*(cj*n+ci) + int(y*n-cj > x*n-ci)
            c = coarse_basis[tid] @ u[coarse_dofs[tid]]
            v = np.array([c[0]+x*c[1]+y*c[2],
                          c[3]+x*c[4]+y*c[5]])
            block += area[k]*w*(grads @ grads.T)
            load += area[k]*w*(grads @ v)
        ii, jj = np.meshgrid(ids, ids, indexing="ij")
        rows.extend(ii.ravel()); cols.extend(jj.ravel()); data.extend(block.ravel())
        rhs[ids] += load
    lap = sparse.coo_matrix((data, (rows, cols)),
                            shape=(nv+len(edges), nv+len(edges))).tocsr()
    phi = np.zeros(len(rhs))
    phi[1:] = spsolve(lap[1:, 1:], rhs[1:])
    return xy, tri, edges, te, phi


def evaluate_neumann_gradient(fine, resolution=600):
    xy, tri, edges, te, phi = fine
    n = int(round(np.sqrt(len(xy)))) - 1
    coords = (np.arange(resolution)+.5)/resolution
    xx, yy = np.meshgrid(coords, coords)
    ii = np.minimum((xx*n).astype(int), n-1)
    jj = np.minimum((yy*n).astype(int), n-1)
    fx, fy = xx*n-ii, yy*n-jj
    tid = 2*(jj*n+ii) + (fy > fx).astype(int)
    gradlam, _ = geometry(xy, tri)
    grad = np.zeros((2, resolution, resolution))
    for k in range(len(tri)):
        mask = tid == k
        if not np.any(mask):
            continue
        pts = np.column_stack((xx[mask], yy[mask]))
        p = xy[tri[k]]
        lam12 = (pts-p[0]) @ np.linalg.inv(np.stack((p[1]-p[0], p[2]-p[0]), axis=-1)).T
        lam = np.column_stack((1-lam12.sum(axis=1), lam12))
        g = gradlam[k]
        values = phi[np.r_[tri[k], len(xy)+te[k]]]
        deriv = np.zeros((len(pts), 2))
        for a in range(3):
            deriv += values[a]*(4*lam[:, a]-1)[:, None]*g[a]
        for index, (a, b) in enumerate(((0, 1), (1, 2), (2, 0)), start=3):
            deriv += 4*values[index]*(lam[:, a, None]*g[b] + lam[:, b, None]*g[a])
        grad[:, mask] = deriv.T
    return grad


def plot_figure(path, title, source, smooth, divergence_free, patch, info,
                smooth_name="P1 smooth", sol_name="Rotgrad P2"):
    plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})
    fig, axes = plt.subplots(2, 4, figsize=(16, 9), dpi=120)
    axes = axes.ravel()
    layers = [source[0], source[1], smooth[0], smooth[1],
              divergence_free[0], divergence_free[1], patch]
    titles = ["BDM field · x", "BDM field · y", f"{smooth_name} · x",
              f"{smooth_name} · y", f"{sol_name} · x",
              f"{sol_name} · y", "Σᵥ (|uᵥ|² + |div uᵥ|²)"]
    scales = [max(np.max(np.abs(layer)) for layer in layers[j:j+2])
              for j in (0, 2, 4)]
    for k, (arr, label) in enumerate(zip(layers, titles)):
        ax = axes[k]
        if k == 6 and patch is None:
            ax.text(.5, .5, "No vertex term", ha="center", va="center",
                    transform=ax.transAxes, fontsize=16)
            ax.set_facecolor("#f3f3f3")
        else:
            vmax = scales[k//2] if k < 6 else None
            im = ax.imshow(arr, origin="lower", extent=(0, 1, 0, 1),
                           cmap="viridis" if k == 6 else "RdBu_r",
                           vmin=0 if k == 6 else -vmax,
                           vmax=None if k == 6 else vmax, interpolation="nearest")
        ax.set_title(label, pad=5)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal")
    axes[7].axis("off")
    axes[7].text(.04, .96, info, ha="left", va="top", transform=axes[7].transAxes,
                 fontsize=13, linespacing=1.6)
    scale_info = (f"Color ranges: source ±{scales[0]:.3g}; "
                  f"smooth ±{scales[1]:.3g};\n"
                  f"divergence-free ±{scales[2]:.3g}")
    if patch is not None:
        scale_info += f"; patch 0–{np.max(patch):.3g}"
    axes[7].text(.04, .08, scale_info,
                 ha="left", va="bottom", transform=axes[7].transAxes,
                 fontsize=11)
    fig.suptitle(title, fontsize=19, y=.986)
    fig.subplots_adjust(left=.018, right=.985, bottom=.035, top=.925,
                        wspace=.065, hspace=.16)
    fig.savefig(path, dpi=120)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "output")
    parser.add_argument("--helmholtz-resolution", type=int, default=40)
    args = parser.parse_args()
    if args.helmholtz_resolution < 10 or args.helmholtz_resolution % 10:
        parser.error("--helmholtz-resolution must be a positive multiple of 10")
    args.output.mkdir(parents=True, exist_ok=True)
    n = 10
    xy, tri, edges, te, _ = mesh(n)
    normal, tri_dofs, basis = bdm_data(xy, tri, edges, te)
    # A basis function attached to a vertex vanishes at the other two
    # vertices of every triangle; this is the requested patch property.
    max_other_vertex = 0.0
    for k, t in enumerate(tri):
        for j, dof in enumerate(tri_dofs[k]):
            owner = edges[dof//2, dof % 2]
            for v in t:
                if v == owner:
                    continue
                x, y = xy[v]
                max_other_vertex = max(max_other_vertex,
                    abs(basis[k, 0, j]+x*basis[k, 1, j]+y*basis[k, 2, j]),
                    abs(basis[k, 3, j]+x*basis[k, 4, j]+y*basis[k, 5, j]))
    assert max_other_vertex < 1e-12, max_other_vertex
    hdiv = bdm_hdiv(xy, tri, tri_dofs, basis)
    h1 = assemble_h1(xy, tri, te, len(xy), 1)
    h2 = assemble_h1(xy, tri, te, len(xy), 2)
    smap, qmap = prolongations(xy, edges, normal)
    pre = preconditioner(smap, qmap, h1, h2, hdiv, edges, len(xy))
    fields, source_data = source_fields(xy, edges, normal, qmap,
                                        return_components=True)
    report = {"mesh": "10x10 square cells, each split into two triangles",
              "bdm_dimension": len(edges)*2, "p1_vector_dimension": 2*len(xy),
              "p2_scalar_dimension": len(xy)+len(edges),
              "max_patch_value_at_other_vertices": max_other_vertex,
              "helmholtz_p2_mesh": args.helmholtz_resolution,
              "noisy_source": {"seed": 20260919,
                               "smooth_h1_norm": source_data["smooth_h1_norm"],
                               "potential_h1_norm": source_data["potential_h1_norm"]},
              "results": {}}
    for name, u in fields.items():
        source, tid, xx, yy = evaluate_bdm_image(u, n, tri_dofs, basis)
        report["results"][name] = {}
        for alpha in (.25, 1., 4., None):
            smooth, sol, local, patches, err, energies = optimal_split(u, pre, alpha)
            smooth_im, _, _, _ = evaluate_bdm_image(smooth, n, tri_dofs, basis)
            sol_im, _, _, _ = evaluate_bdm_image(sol, n, tri_dofs, basis)
            patch = patch_energy_image(patches, tri_dofs, basis, tid, xx, yy) \
                    if patches is not None else None
            key = "no_local" if alpha is None else f"alpha_{alpha:g}"
            label = "P1 + rotgrad P2, no local term" if alpha is None else f"HX optimal split · α = {alpha:g}"
            info = ("Unit square · 10×10 cells\n200 triangles · BDM1\n"
                    f"Reconstruction error: {err:.2e}\n"
                    f"Minimum product norm²: {energies['weighted_total']:.4g}\n"
                    f"P1 H¹ norm²: {energies['smooth_h1_squared']:.4g}\n"
                    f"P2 H¹ norm²: {energies['potential_h1_squared']:.4g}\n"
                    f"Vertex H(div) norm²: {energies['local_hdiv_squared']:.4g}")
            plot_figure(args.output / f"{name}_{key}.png", f"{name.replace('_', ' ').title()} · {label}",
                        source, smooth_im, sol_im, patch, info)
            report["results"][name][key] = {"reconstruction_l2_dof": err,
                                              **energies}
            print(name, key, f"error={err:.2e}",
                  f"energy={energies['weighted_total']:.6g}", flush=True)
        fine = neumann_helmholtz(u, n, tri_dofs, basis,
                                 args.helmholtz_resolution)
        gradient = evaluate_neumann_gradient(fine)
        divergence_free = source - gradient
        info = ("Continuous Neumann Helmholtz split\n"
                f"P2 Galerkin mesh: {args.helmholtz_resolution}×{args.helmholtz_resolution}\n"
                "Curl-free part: ∇φ\n"
                "∫ ∇φ·∇ψ = ∫ u·∇ψ\n"
                "Remainder: u − ∇φ\n"
                "Remainder has zero normal flux\nin the weak limit")
        plot_figure(args.output / f"{name}_helmholtz.png",
                    f"{name.replace('_', ' ').title()} · Neumann Helmholtz split",
                    source, gradient, divergence_free, None, info,
                    smooth_name="Curl-free ∇φ", sol_name="Divergence-free")
        report["results"][name]["helmholtz"] = {"method": "P2 Neumann Galerkin approximation"}
        print(name, "helmholtz", flush=True)
    (args.output / "metrics.json").write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    main()
