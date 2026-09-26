#!/usr/bin/env python3
"""Full product-energy constants against the H(div) norm, alpha = 1.

CoptHX: lambda_max(B^{-1}, A) with B = bs + bq + bl, the exact discrete HX optimum.
Cavg:   one quadratic form ||J_h s_r||^2 + ||Pi_Q q_r||^2 + sum_patch ||local||^2,
        with (s_r, q_r) the continuous Galerkin minimizer at spline refinement r.
Ccont:  ||s_r||^2 + ||q_r||^2 of the continuous minimizer (no transfer).

The spline correction E = E' B0^{-1} is never formed: E' = -G^{-1} K is solved in
column chunks and only its small images are kept.
"""
import argparse
import json
import resource
import time

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import splu

import averaged_transfer
import maximize_patch_ratio as search
import plot_decompositions as base
import sparse_continuous
from paths import AGENT_ARTIFACTS


def eig_end(a, b, top=True):
    k = len(b)-1 if top else 0
    val, vec = linalg.eigh((a+a.T)/2, b, subset_by_index=(k, k), check_finite=False)
    u = vec[:, 0]
    return float(val[0]), (u if u[np.argmax(np.abs(u))] > 0 else -u)


def assemble_loads(xy, tri, te, potential, refinement):
    """Sparse copy of sparse_continuous.assemble loads (same formulas)."""
    nv, ne, m = len(xy), int(te.max())+1, len(potential.modes)
    gradlam, area = base.geometry(xy, tri)
    lam, weights = averaged_transfer.subtriangle_rule(refinement)
    blocks_s, blocks_q = [], []
    for k, t in enumerate(tri):
        val, gx, gy, hxx, hxy, hyy = sparse_continuous.evaluate(potential, lam @ xy[t])
        w = 2*area[k]*weights
        v, x, y, xx, cross, yy = (mat.T.multiply(w).tocsr() for mat in (val, gx, gy, hxx, hxy, hyy))
        g = gradlam[k]
        g0, g1 = np.broadcast_to(g[:, 0], lam.shape), np.broadcast_to(g[:, 1], lam.shape)
        c0 = y @ lam + cross @ g0 + yy @ g1
        c1 = -x @ lam - xx @ g0 - cross @ g1
        shape = np.empty((len(lam), 6)); grad = np.empty((len(lam), 6, 2))
        shape[:, :3] = lam*(2*lam-1)
        grad[:, :3] = (4*lam-1)[:, :, None]*g
        for j, (a, b) in enumerate(((0, 1), (1, 2), (2, 0)), 3):
            shape[:, j] = 4*lam[:, a]*lam[:, b]
            grad[:, j] = 4*(lam[:, a, None]*g[b]+lam[:, b, None]*g[a])
        cqk = -v @ shape - x @ grad[:, :, 0] - y @ grad[:, :, 1]
        rows = np.unique(val.indices)
        for mat, cols, out in ((c0, 2*t, blocks_s), (c1, 2*t+1, blocks_s),
                               (cqk, np.r_[t, nv+te[k]], blocks_q)):
            sub = np.asarray(mat)[rows]
            out.append((np.repeat(rows, len(cols)), np.tile(cols, len(rows)), sub.ravel()))
    def build(blocks, ncol):
        i, j, v = (np.concatenate(z) for z in zip(*blocks))
        return sparse.csr_matrix((v, (i, j)), shape=(m, ncol))
    return build(blocks_s, 2*nv), build(blocks_q, nv+ne)


def setup(n):
    xy, tri, edges, te, _ = base.mesh(n)
    normal, td, basis = base.bdm_data(xy, tri, edges, te)
    hdiv = base.bdm_hdiv(xy, tri, td, basis)
    h1 = base.assemble_h1(xy, tri, te, len(xy), 1)
    h2 = base.assemble_h1(xy, tri, te, len(xy), 2)
    smap, qmap = base.prolongations(xy, edges, normal)
    pre = base.preconditioner(smap, qmap, h1, h2, hdiv, edges, len(xy))
    a = pre['av']
    w = np.zeros_like(a)
    for ids, _ in pre['fact']:
        w[np.ix_(ids, ids)] = a[np.ix_(ids, ids)]
    eye = np.eye(len(a))
    b0inv = linalg.solve(pre['bs']+pre['bq'], eye, assume_a='sym', check_finite=False)
    binv = linalg.solve(pre['bs']+pre['bq']+pre['bl'], eye, assume_a='sym', check_finite=False)
    ss = linalg.cho_solve(pre['sf'], pre['sd'].T, check_finite=False)
    qq = linalg.cho_solve(pre['qf'], pre['qd'].T, check_finite=False)
    return dict(n=n, xy=xy, tri=tri, edges=edges, te=te, pre=pre, a=a, w=w,
                b0inv=(b0inv+b0inv.T)/2, binv=(binv+binv.T)/2, ss=ss, qq=qq,
                hdiv_sparse_diff=float(abs(hdiv-sparse.csr_matrix(a)).max()))


def copt_hx(st):
    pre, a, binv = st['pre'], st['a'], st['binv']
    val, u = eig_end(binv, a)
    b = pre['bs']+pre['bq']+pre['bl']
    split = base.optimal_split(u, pre, 1.0)
    energies = split[-1] if isinstance(split[-1], dict) else None
    lagrange = linalg.solve(b, u, assume_a='sym')
    s = linalg.cho_solve(pre['sf'], pre['sd'].T @ lagrange)
    q = linalg.cho_solve(pre['qf'], pre['qd'].T @ lagrange)
    local = np.zeros_like(u); le = 0.0
    for ids, fac in pre['fact']:
        c = linalg.cho_solve(fac, lagrange[ids]); local[ids] = c
        le += float(c @ a[np.ix_(ids, ids)] @ c)
    ua = float(u @ a @ u)
    es, eq = float(s @ pre['hs'] @ s), float(q @ pre['hq'] @ q)
    ev = binv @ u - val*(a @ u)
    rec = dict(CoptHX=val, rayleigh=float(u @ binv @ u)/ua,
               split_energy_ratio=(es+eq+le)/ua,
               split_components=dict(s_h1=es/ua, q_h1=eq/ua, patches=le/ua),
               reconstruction_relative=float(np.linalg.norm(pre['sd'] @ s+pre['qd'] @ q+local-u)/np.linalg.norm(u)),
               eigen_residual_relative=float(np.linalg.norm(ev)/np.linalg.norm(binv @ u)),
               binv_symmetry=float(abs(binv-binv.T).max()/abs(binv).max()))
    if energies is not None:
        rec['optimal_split_energies'] = {k: v/ua for k, v in energies.items() if isinstance(v, float)}
    # No-local optimum for reference: lambda_max(B0^{-1}, A).
    rec['C_nolocal_discrete'] = eig_end(st['b0inv'], a)[0]
    return rec, u


def continuous_maps(st, r, chunk):
    n, xy, tri, te, edges, pre = st['n'], st['xy'], st['tri'], st['te'], st['edges'], st['pre']
    potential = search.TensorSplinePotential(n*r)
    t0 = time.monotonic()
    gram = sparse_continuous.gram_matrix(potential)
    gq = sparse_continuous.gram_matrix(potential, potential_only=True)
    cs, cq = assemble_loads(xy, tri, te, potential, r)
    nodal, scalar = sparse_continuous.moments(potential, xy, edges, r)
    jmap = averaged_transfer.averaged_rotgrad(potential, xy, tri, r)
    t1 = time.monotonic()
    lu = splu(gram.tocsc())
    t2 = time.monotonic()
    N = len(st['a'])
    ss, qq = st['ss'], st['qq']
    out = dict(jE=np.zeros((jmap.shape[0], N)), sE=np.zeros((scalar.shape[0], N)),
               nE=np.zeros((nodal.shape[0], N)), KE=np.zeros((N, N)), c1E=np.zeros(N))
    cq1 = cq @ np.ones(cq.shape[1])
    res = ident = 0.0
    for c0 in range(0, N, chunk):
        sl = slice(c0, min(N, c0+chunk))
        k = cs @ ss[:, sl] + cq @ qq[:, sl]
        e = -lu.solve(k)
        ge = gram @ e
        res = max(res, np.linalg.norm(ge+k)/np.linalg.norm(k))
        out['jE'][:, sl] = jmap @ e
        out['sE'][:, sl] = scalar @ e
        out['nE'][:, sl] = nodal @ e
        out['c1E'][sl] = cq1 @ e
        # K^T E' (N x chunk); E'^T G E' = -K^T E' is checked on the diagonal block.
        out['KE'][:, sl] = ss.T @ (cs.T @ e) + qq.T @ (cq.T @ e)
        ident = max(ident, np.abs(e.T @ ge + out['KE'][sl, sl]).max()/np.abs(out['KE'][sl, sl]).max())
        del k, e, ge
    t3 = time.monotonic()
    out.update(spline_dimension=gram.shape[0], stationarity_residual=float(res),
               energy_identity_block_error=float(ident),
               # lu.nnz is SuperLU storage; lu.L/lu.U would copy the factor.
               lu_storage_nnz=int(lu.nnz),
               timing=dict(assemble=t1-t0, factor=t2-t1, solve=t3-t2))
    ops = dict(lu=lu, gram=gram, gq=gq, cs=cs, cq=cq, jmap=jmap, scalar=scalar)
    return out, ops


def psi_of(st, ops, u):
    """Spline coefficients of the continuous correction for one input u."""
    lam0 = st['b0inv'] @ u
    k = ops['cs'] @ (st['ss'] @ lam0) + ops['cq'] @ (st['qq'] @ lam0)
    return -ops['lu'].solve(k), st['ss'] @ lam0, st['qq'] @ lam0


def analyse(st, out, ops, r, previous):
    pre, a, w, b0inv, binv = st['pre'], st['a'], st['w'], st['b0inv'], st['binv']
    sd, qd, hs, hq = pre['sd'], pre['qd'], pre['hs'], pre['hq']
    s0, q0 = st['ss'] @ b0inv, st['qq'] @ b0inv
    je, se, ne_ = out['jE'] @ b0inv, out['sE'] @ b0inv, out['nE'] @ b0inv
    ms, mq, lmap = s0+je, q0-se, qd @ se - sd @ je
    fs, fq, fl = ms.T @ hs @ ms, mq.T @ hq @ mq, lmap.T @ w @ lmap
    f = fs+fq+fl
    cavg, u = eig_end(f, a)
    ua = float(u @ a @ u)
    ones = np.ones(hq.shape[0]); mvec = hq @ ones; area = float(ones @ mvec)
    mq0 = mq - np.outer(ones, mvec @ mq)/area
    fq0 = mq0.T @ hq @ mq0
    cavg0, u0 = eig_end(fs+fq0+fl, a)
    lnod = qd @ se - sd @ ne_
    # Continuous minimizer form and the two direct spline evaluations first, so the
    # spline factorization can be released before the remaining dense eigensolves.
    c = b0inv + b0inv @ out['KE'] @ b0inv
    c = (c+c.T)/2
    ccont, uc = eig_end(c, a)
    meanrow = mvec @ q0 + out['c1E'] @ b0inv
    def direct(v):
        psi, sc, qc = psi_of(st, ops, v)
        gp, gqp = ops['gram'] @ psi, ops['gq'] @ psi
        e_s = sc @ hs @ sc + 2*psi @ (ops['cs'] @ sc) + psi @ (gp-gqp)
        e_q = qc @ hq @ qc + 2*psi @ (ops['cq'] @ qc) + psi @ gqp
        va = float(v @ a @ v)
        return dict(s_h1=float(e_s)/va, q_h1=float(e_q)/va, form=float(v @ c @ v)/va,
                    q_mean=float(meanrow @ v)/np.sqrt(va),
                    J_psi_error=float(np.abs(ops['jmap'] @ psi - je @ v).max()))
    direct_uc, direct_u = direct(uc), direct(u)
    ops.clear()
    rec = dict(refinement=r, spline_dimension=out['spline_dimension'],
               stationarity_residual=out['stationarity_residual'],
               energy_identity_block_error=out['energy_identity_block_error'],
               lu_storage_nnz=out['lu_storage_nnz'], timing=out['timing'])
    rec['Cavg'] = cavg
    rec['Cavg_rayleigh'] = float(u @ f @ u)/ua
    rec['Cavg_components_at_maximizer'] = dict(
        s_h1=float(u @ fs @ u)/ua, q_h1=float(u @ fq @ u)/ua, patches=float(u @ fl @ u)/ua)
    rec['Cavg_eigen_residual_relative'] = float(np.linalg.norm(f @ u - cavg*(a @ u))/np.linalg.norm(f @ u))
    rec['Cavg_mean_subtracted'] = cavg0
    rec['Cavg_mean_subtracted_components_at_maximizer'] = dict(
        s_h1=float(u0 @ fs @ u0)/float(u0 @ a @ u0), q_h1=float(u0 @ fq0 @ u0)/float(u0 @ a @ u0),
        patches=float(u0 @ fl @ u0)/float(u0 @ a @ u0))
    rec['PiQ_mean_operator_norm'] = float(np.sqrt(max(0, eig_end(np.outer(mvec @ mq, mvec @ mq)/area, a)[0])))
    rec['rotgrad_of_constant_norm'] = float(np.abs(qd @ ones).max())
    rec['domain_area'] = area
    rec['smooth_part_operator_norm_squared'] = eig_end(fs, a)[0]
    rec['potential_part_operator_norm_squared'] = eig_end(fq, a)[0]
    rec['patch_only_averaged'] = eig_end(fl, a)[0]
    rec['patch_only_nodal'] = eig_end(lnod.T @ w @ lnod, a)[0]
    rec['reconstruction_max_error'] = float(np.abs(sd @ ms + qd @ mq + lmap - np.eye(len(a))).max())
    rec['dominance_min_eig_Favg_minus_Binv'] = eig_end(f-binv, a, top=False)[0]
    rec['Ccont'] = ccont
    rec['Ccont_eigen_residual_relative'] = float(np.linalg.norm(c @ uc - ccont*(a @ uc))/np.linalg.norm(c @ uc))
    rec['nolocal_minus_Ccont_min_eig'] = eig_end(b0inv-c, a, top=False)[0]
    rec['continuous_q_mean_operator_norm'] = float(np.sqrt(max(0, eig_end(np.outer(meanrow, meanrow), a)[0])))
    rec['Ccont_mean_subtracted'] = eig_end(c - np.outer(meanrow, meanrow)/area, a)[0]
    rec['Ccont_direct_at_maximizer'] = direct_uc
    rec['continuous_energies_at_Cavg_maximizer'] = direct_u
    maps = dict(ms=ms, mq=mq, l=lmap, c=c, u=u, uc=uc)
    if previous is not None:
        dms, dmq, dl = ms-previous['ms'], mq-previous['mq'], lmap-previous['l']
        ch = lambda x: float(np.sqrt(max(0, eig_end(x, a)[0])))
        rec['compared_to_refinement'] = previous['r']
        rec['product_map_change_norm'] = ch(dms.T @ hs @ dms + dmq.T @ hq @ dmq + dl.T @ w @ dl)
        rec['map_change_components'] = dict(smooth=ch(dms.T @ hs @ dms), potential=ch(dmq.T @ hq @ dmq),
                                            patches=ch(dl.T @ w @ dl))
        pu = previous['u']
        rec['previous_Cavg_maximizer_ratio_now'] = float(pu @ f @ pu)/float(pu @ a @ pu)
        dc = previous['c'] - c
        rec['Ccont_form_decrease_max'] = eig_end(dc, a)[0]
        rec['Ccont_form_decrease_min'] = eig_end(dc, a, top=False)[0]
        rec['continuous_map_change_norm'] = float(np.sqrt(max(0, rec['Ccont_form_decrease_max'])))
    return rec, maps


def peak_rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mesh', type=int, default=5)
    parser.add_argument('--refinements', nargs='*', type=int, default=[8, 16, 32])
    parser.add_argument('--chunk', type=int, default=256)
    args = parser.parse_args()
    path = AGENT_ARTIFACTS/f'full_constants_n{args.mesh}.json'
    data = json.loads(path.read_text()) if path.exists() else dict(mesh=args.mesh, refinements=[])
    t = time.monotonic()
    st = setup(args.mesh)
    data['bdm_dimension'] = len(st['a'])
    data['hdiv_dense_consistency'] = st['hdiv_sparse_diff']
    copt, u = copt_hx(st)
    copt['elapsed_seconds'] = time.monotonic()-t
    copt['peak_rss_gb'] = peak_rss_gb()
    data['CoptHX'] = copt
    np.save(AGENT_ARTIFACTS/f'full_constants_copt_maximizer_n{args.mesh}.npy', u)
    path.write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps(copt), flush=True)
    previous = None
    for r in args.refinements:
        t = time.monotonic()
        out, ops = continuous_maps(st, r, args.chunk)
        rec, maps = analyse(st, out, ops, r, previous)
        rec['CoptHX_le_Cavg'] = bool(copt['CoptHX'] <= rec['Cavg']*(1+1e-12))
        rec['elapsed_seconds'] = time.monotonic()-t
        rec['peak_rss_gb'] = peak_rss_gb()
        np.save(AGENT_ARTIFACTS/f'full_constants_avg_maximizer_n{args.mesh}_ref{r}.npy', maps['u'])
        np.save(AGENT_ARTIFACTS/f'full_constants_cont_maximizer_n{args.mesh}_ref{r}.npy', maps['uc'])
        data['refinements'] = sorted([x for x in data['refinements'] if x['refinement'] != r]+[rec],
                                     key=lambda x: x['refinement'])
        path.write_text(json.dumps(data, indent=2)+'\n')
        print(json.dumps(rec), flush=True)
        previous = dict(maps, r=r)
        del out, ops


if __name__ == '__main__':
    main()
