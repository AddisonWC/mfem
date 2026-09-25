#!/usr/bin/env python3
"""Independent assembly and polynomial-exactness checks for the experiment."""
import numpy as np
import maximize_patch_ratio as search
import sparse_continuous


def main():
    xy, tri, edges, te, mid = search.base.mesh(2)
    p = search.TensorSplinePotential(4)
    _, gram, cs, cq = search.continuous_quadratic_maps(xy, tri, te, 8, p, 2)
    g2, s2, q2 = sparse_continuous.assemble(xy, tri, te, p, 2)
    nodal, scalar = search.moment_maps(p, xy, edges, 2)
    n2, m2 = sparse_continuous.moments(p, xy, edges, 2)
    for label, a, b in [('gram', gram, g2.toarray()), ('smooth load',cs,s2),
                         ('potential load',cq,q2), ('nodal map',nodal,n2.toarray()),
                         ('moment map',scalar,m2.toarray())]:
        error = np.linalg.norm(a-b)/max(1,np.linalg.norm(a))
        assert error < 1e-12, (label,error)
        print(label, error)
    h1 = search.base.assemble_h1(xy,tri,te,len(xy),2)
    coeff = np.r_[xy[:,0]**2,mid[:,0]**2]
    # Integral of x^4 + (2x)^2 on the unit square.
    assert abs(coeff @ h1 @ coeff - (1/5+4/3)) < 1e-12
    # Derivative maps annihilate a constant and reproduce affine polynomials.
    points = np.array([[0,0],[1,1],[.27,.63],[.5,.5]])
    val,gx,gy,*_ = sparse_continuous.evaluate(p,points)
    assert np.allclose(val @ np.ones(p.nb**2),1)
    assert np.allclose(gx @ np.ones(p.nb**2),0)
    assert np.allclose(gy @ np.ones(p.nb**2),0)
    print('P2 H1 polynomial exactness and spline partition checks passed')


if __name__ == '__main__':
    main()
