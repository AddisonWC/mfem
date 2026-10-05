#!/usr/bin/env python3
"""Full BDM1 HX spectra on domains with holes, for div-div + tau mass.

The domain is the unit square cut into an n x n grid of diagonally split
squares, with square blocks of cells removed. Each hole adds one harmonic
H(div) field (flux around the hole), so the cohomology has dimension equal to
the number of holes.

The additive preconditioner is HX (7.3) with exact subspace solves:

    B = R + P (L + tau M)^{-1} P^T + Q (Q^T A Q)^{-1} Q^T,

where P interpolates continuous vector P1, Q = rot from continuous P2 (one
gauge node removed), and R is Jacobi or exact vertex-patch solves. L is the
unweighted vector H1 seminorm. The "galerkin" auxiliary replaces
(L + tau M)^{-1} by (P^T A P)^{-1}, as in the matrix form used by ADS.
"""

import argparse
import json

import numpy as np
from scipy import linalg, sparse

from paths import AGENT_ARTIFACTS
from plot_decompositions import (QPTS, QWTS, assemble_h1, bdm_data, geometry,
                                 prolongations)


def holed_mesh(n, holes):
    """Grid mesh with the cells [i0, i1) x [j0, j1) of each hole removed."""
    removed = lambda i, j: any(i0 <= i < i1 and j0 <= j < j1
                               for i0, i1, j0, j1 in holes)
    used = {}
    xy = []
    def vertex(i, j):
        if (i, j) not in used:
            used[(i, j)] = len(xy)
            xy.append((i / n, j / n))
        return used[(i, j)]
    tri = []
    for j in range(n):
        for i in range(n):
            if removed(i, j):
                continue
            a, b, c, d = (vertex(i, j), vertex(i + 1, j),
                          vertex(i + 1, j + 1), vertex(i, j + 1))
            tri.extend(((a, b, c), (a, c, d)))
    xy = np.asarray(xy, dtype=float)
    tri = np.asarray(tri, dtype=int)
    edge_ids, edges, te = {}, [], []
    for t in tri:
        local = []
        for x, y in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            pair = tuple(sorted((x, y)))
            if pair not in edge_ids:
                edge_ids[pair] = len(edges)
                edges.append(pair)
            local.append(edge_ids[pair])
        te.append(local)
    return xy, tri, np.asarray(edges, dtype=int), np.asarray(te, dtype=int)


def bdm_div_mass(xy, tri, tri_dofs, basis):
    """Separate BDM1 div-div and mass matrices (3-point rule is exact)."""
    _, area = geometry(xy, tri)
    nd = int(tri_dofs.max()) + 1
    rows, cols, dd, mm = [], [], [], []
    for k, t in enumerate(tri):
        mass = np.zeros((6, 6))
        div = basis[k, 1] + basis[k, 5]
        for w, lam in zip(QWTS, QPTS):
            x, y = lam @ xy[t]
            v = np.stack((basis[k, 0] + x*basis[k, 1] + y*basis[k, 2],
                          basis[k, 3] + x*basis[k, 4] + y*basis[k, 5]))
            mass += area[k] * w * (v.T @ v)
        ids = tri_dofs[k]
        ii, jj = np.meshgrid(ids, ids, indexing="ij")
        rows.extend(ii.ravel()); cols.extend(jj.ravel())
        dd.extend((area[k] * np.outer(div, div)).ravel())
        mm.extend(mass.ravel())
    shape = (nd, nd)
    return (sparse.coo_matrix((dd, (rows, cols)), shape=shape).toarray(),
            sparse.coo_matrix((mm, (rows, cols)), shape=shape).toarray())


def p1_mass_stiffness(xy, tri, te):
    """Split the P1 H1 matrix from assemble_h1 into mass and stiffness."""
    gradlam, area = geometry(xy, tri)
    nv = len(xy)
    rows, cols, kk = [], [], []
    for k, t in enumerate(tri):
        ii, jj = np.meshgrid(t, t, indexing="ij")
        rows.extend(ii.ravel()); cols.extend(jj.ravel())
        kk.extend((area[k] * gradlam[k] @ gradlam[k].T).ravel())
    stiff = sparse.coo_matrix((kk, (rows, cols)), shape=(nv, nv)).toarray()
    mass = assemble_h1(xy, tri, te, nv, 1).toarray() - stiff
    return stiff, mass


def harmonic_dimension(div, mass, qmap):
    """dim ker(div) - dim range(rot), computed from the assembled matrices."""
    kernel = int(np.sum(linalg.eigvalsh(div, mass) < 1e-9))
    return kernel - np.linalg.matrix_rank(qmap)


def spectrum(div, mass, smap, qmap, stiff1, mass1, edges, tau, smoother,
             aux):
    a = div + tau * mass
    if smoother == "jacobi":
        b = np.diag(1 / np.diag(a))
    else:
        b = np.zeros_like(a)
        patches = [[] for _ in range(int(edges.max()) + 1)]
        for e, pair in enumerate(edges):
            for side, v in enumerate(pair):
                patches[v].append(2*e + side)
        for ids in patches:
            ids = np.asarray(ids, dtype=int)
            b[np.ix_(ids, ids)] = linalg.inv(a[np.ix_(ids, ids)])
    if aux == "hx":
        h1 = np.kron(stiff1 + tau * mass1, np.eye(2))
    else:
        h1 = smap.T @ a @ smap
    b += smap @ linalg.solve(h1, smap.T, assume_a="pos")
    qr = qmap[:, 1:]  # rot annihilates constants; drop one gauge node.
    b += qr @ linalg.solve(qr.T @ a @ qr, qr.T, assume_a="pos")
    # With A = L L^T, BA is similar to the symmetric L^T B L.
    chol = linalg.cholesky(a, lower=True)
    return linalg.eigvalsh(chol.T @ b @ chol)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-n", type=int, default=16)
    parser.add_argument("--holes", default="1",
                        help="0, 1, or 2 square holes")
    parser.add_argument("--smoother", choices=("jacobi", "patch"),
                        default="patch")
    parser.add_argument("--aux", choices=("hx", "galerkin"), default="hx")
    parser.add_argument("--taus", default="1e-6,1e-4,1e-2,1,1e2")
    parser.add_argument("--show", type=int, default=5,
                        help="number of smallest and largest values to print")
    args = parser.parse_args()

    n = args.n
    q = n // 8
    layouts = {"0": [], "1": [(3*q, 5*q, 3*q, 5*q)],
               "2": [(q, 3*q, 3*q, 5*q), (5*q, 7*q, 3*q, 5*q)]}
    xy, tri, edges, te = holed_mesh(n, layouts[args.holes])
    normal, tri_dofs, basis = bdm_data(xy, tri, edges, te)
    div, mass = bdm_div_mass(xy, tri, tri_dofs, basis)
    smap, qmap = (m.toarray() for m in prolongations(xy, edges, normal))
    stiff1, mass1 = p1_mass_stiffness(xy, tri, te)
    dim_h = harmonic_dimension(div, mass, qmap)
    print(f"# n={n} holes={args.holes} bdm_dofs={len(div)} "
          f"harmonic_dim={dim_h} smoother={args.smoother} aux={args.aux}")
    records = []
    for tau in (float(t) for t in args.taus.split(",")):
        lam = spectrum(div, mass, smap, qmap, stiff1, mass1, edges, tau,
                       args.smoother, args.aux)
        k = args.show
        print(f"tau={tau:8.1e} kappa={lam[-1]/lam[0]:10.3e} "
              f"kappa_drop{dim_h}={lam[-1]/lam[dim_h]:8.3f}  "
              f"smallest={np.array2string(lam[:k], precision=4)} "
              f"largest={np.array2string(lam[-k:], precision=4)}")
        records.append(dict(tau=tau, smallest=lam[:k].tolist(),
                            largest=lam[-k:].tolist()))
    out = AGENT_ARTIFACTS / (f"cohomology_n{n}_holes{args.holes}_"
                             f"{args.smoother}_{args.aux}.json")
    out.write_text(json.dumps(dict(n=n, holes=args.holes,
                                   harmonic_dim=int(dim_h),
                                   smoother=args.smoother, aux=args.aux,
                                   spectra=records), indent=1))


if __name__ == "__main__":
    main()
