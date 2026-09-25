"""Sparse, exactly integrated H2 tensor-spline correction operators."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.interpolate import BSpline

import plot_decompositions as base


def evaluate(potential, points):
    """Six sparse value/derivative matrices, at most 16 entries per row."""
    pts = np.asarray(points).reshape(-1, 2)
    nb = potential.nb
    # All derivatives have the same four possibly nonzero basis functions.
    ix = np.minimum((pts[:, 0]*potential.n).astype(int), potential.n-1)
    iy = np.minimum((pts[:, 1]*potential.n).astype(int), potential.n-1)
    ix = ix[:, None] + np.arange(4)
    iy = iy[:, None] + np.arange(4)
    rows = np.arange(len(pts))[:, None]
    def univariate(x, indices):
        knots = potential.knots
        transform = sparse.eye(nb, format='csr')
        values = []
        for k in range(3):
            matrix = sparse.csr_matrix(BSpline.design_matrix(x, knots, 3-k) @ transform)
            selected = matrix[rows, indices]
            values.append(selected.toarray() if sparse.issparse(selected) else np.asarray(selected))
            if k < 2:
                size = len(knots)-(3-k)-1
                scale = (3-k)/(knots[4-k:4-k+size-1]-knots[1:size])
                derivative = sparse.diags((-scale,scale),(0,1),shape=(size-1,size))
                transform = derivative @ transform
                knots = knots[1:-1]
        return values
    bx, by = univariate(pts[:,0],ix), univariate(pts[:,1],iy)
    cols = (ix[:, :, None]*nb+iy[:, None, :]).reshape(-1)
    rows = np.repeat(np.arange(len(pts)), 16)
    return tuple(sparse.csr_matrix(((bx[a][:, :, None]*by[b][:, None, :]).ravel(),
                                    (rows, cols)), shape=(len(pts), nb*nb))
                 for a, b in ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)))


def gram_matrix(potential, potential_only=False):
    z, w = leggauss(4)
    x = ((np.arange(potential.n)[:, None]+(z+1)/2)/potential.n).ravel()
    w = np.tile(w/(2*potential.n), potential.n)
    mats = []
    for k in range(3):
        b = sparse.csr_matrix(potential.spline(x, nu=k))
        mats.append(b.T @ b.multiply(w[:, None]))
    m, d, h = mats
    if potential_only:
        return (sparse.kron(m,m)+sparse.kron(d,m)+sparse.kron(m,d)).tocsc()
    return (sparse.kron(m, m) + 2*sparse.kron(d, m) + 2*sparse.kron(m, d)
            + sparse.kron(h, m) + 2*sparse.kron(d, d) + sparse.kron(m, h)).tocsc()


def assemble(xy, tri, te, potential, refinement):
    """Integrate loads on subtriangles aligned with every spline knot."""
    nv, ne, m = len(xy), int(te.max())+1, len(potential.modes)
    cs, cq = np.zeros((m, 2*nv)), np.zeros((m, nv+ne))
    gradlam, area = base.geometry(xy, tri)
    z, wz = leggauss(5)
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
            a0, b0, c0 = bary(i,j), bary(i+1,j), bary(i,j+1)
            subs.append(np.stack((a0,b0,c0)))
            if i+j < r-1:
                subs.append(np.stack((b0,bary(i+1,j+1),c0)))
    lam = np.concatenate([ref @ sub for sub in subs])
    weights = np.tile(weights, r*r)/(r*r)
    for k, t in enumerate(tri):
        val, gx, gy, hxx, hxy, hyy = evaluate(potential, lam @ xy[t])
        w = 2*area[k]*weights
        weighted = lambda mat: mat.T.multiply(w)
        v, x, y, xx, cross, yy = map(weighted, (val,gx,gy,hxx,hxy,hyy))
        g = gradlam[k]
        cs[:, 2*t] += y @ lam + cross @ np.broadcast_to(g[:,0], lam.shape) + yy @ np.broadcast_to(g[:,1], lam.shape)
        cs[:, 2*t+1] += -x @ lam - xx @ np.broadcast_to(g[:,0], lam.shape) - cross @ np.broadcast_to(g[:,1], lam.shape)
        shape = np.empty((len(lam),6)); grad = np.empty((len(lam),6,2))
        shape[:,:3] = lam*(2*lam-1)
        grad[:,:3] = (4*lam-1)[:,:,None]*g
        for j, (a,b) in enumerate(((0,1),(1,2),(2,0)),3):
            shape[:,j] = 4*lam[:,a]*lam[:,b]
            grad[:,j] = 4*(lam[:,a,None]*g[b]+lam[:,b,None]*g[a])
        cq[:, np.r_[t,nv+te[k]]] += -v @ shape - x @ grad[:,:,0] - y @ grad[:,:,1]
    return gram_matrix(potential), cs, cq


def moments(potential, xy, edges, refinement):
    val, gx, gy, *_ = evaluate(potential, xy)
    nodal = sparse.vstack([gy, -gx]).tocsr()[np.column_stack((np.arange(len(xy)),len(xy)+np.arange(len(xy)))).ravel()]
    z, w = leggauss(4)
    t = ((np.arange(refinement)[:,None]+(z+1)/2)/refinement).ravel()
    w = np.tile(w/(2*refinement),refinement)
    pts = xy[edges[:,0],None,:]+t[None,:,None]*(xy[edges[:,1]]-xy[edges[:,0]])[:,None,:]
    ev = evaluate(potential, pts)[0]
    averaging = sparse.kron(sparse.eye(len(edges)), sparse.csr_matrix(w[None,:]))
    mid = 1.5*(averaging @ ev)-.25*(val[edges[:,0]]+val[edges[:,1]])
    return nodal, sparse.vstack([val,mid]).tocsr()
