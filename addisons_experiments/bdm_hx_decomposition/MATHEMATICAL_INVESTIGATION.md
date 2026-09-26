# Exact-minimizer stability and alternatives to nodal interpolation

September 25, 2026. This note concerns the unit square, its uniform square
grids split along a fixed diagonal, triangular BDM1, and the full H1 norms.
No essential boundary conditions are imposed. The conclusions below distinguish
the exact continuous minimizer from its finite-dimensional approximations.

## Main conclusion

For the **exact** minimizer

\[
 (s,q)=\arg\min_{u_h=s+R q}
       \bigl(\|s\|_{H^1}^2+\|q\|_{H^1}^2\bigr),
 \qquad Rq=(\partial_yq,-\partial_xq),
\]

the original construction with canonical BDM edge moments and nodal P1
interpolation has a mesh-independent bound

\[
 \sum_v\|u_v\|_{H(\mathrm{div})}^2
 \le C\|u_h\|_{H(\mathrm{div})}^2.
 \tag{1}
\]

Thus the requested supremum cannot blow up as these meshes are refined.
This is a qualitative bound, not a determination of the best constant or
proof that the previously observed values near 0.37 are converged.
The proof uses the PDE satisfied by the exact minimizer, not boundedness
of nodal interpolation on arbitrary H1 functions (which is false).

Replacing nodal interpolation by an H1-stable local averaging projection
gives a simpler proof and makes the patch map continuous in the norm in
which the continuous optimization converges. Neither the BDM moment
interpolant nor the potential transfer needs to be replaced for that result.

## 1. Existence, normalization, and the hidden Stokes problem

Write Z for the closed subspace of divergence-free vectors in L2. On the
simply connected square there is a bounded stream-function map

\[
 T:Z\longrightarrow H^1(\Omega)\cap L^2_0(\Omega),
 \qquad R T w=w,\qquad \|Tw\|_{H^1}\le C\|w\|_{L^2}.
\]

Let g=div u_h. A Dirichlet Poisson solution Delta phi=g on the square
provides an H1 vector lifting s_bar=grad phi with
\(\|s_{\rm bar}\|_{H^1}\le C\|g\|_{L^2}\). Taking
q_bar=T(u_h-s_bar) establishes a uniformly bounded admissible pair.
The feasible set is closed and the objective is strictly convex and
coercive, so the minimizer exists, is unique, and satisfies

\[
 \int_\Omega q=0,\qquad
 \|s\|_{H^1}+\|q\|_{H^1}\le C\|u_h\|_{H(\mathrm{div})}.
 \tag{2}
\]

Let T* denote the L2 adjoint of T viewed as a map Z -> L2_0, with values
in Z. For a divergence-free variation v in H1, the corresponding potential
variation is -Tv. Since R is an isometry from gradients to vectors,
stationarity gives

\[
 (\nabla s,\nabla v)=(F,v),\qquad
 F=u_h-2s+T^*q\in L^2.
\]

The divergence map H1^2 -> L2 has a bounded right inverse on the square.
Consequently there is p in L2 such that, for **every** v in H1^2,

\[
 (\nabla s,\nabla v)-(p,\mathrm{div}\,v)=(F,v),
 \qquad\mathrm{div}\,s=g.
 \tag{3}
\]

This is a Stokes problem with full-gradient natural traction
\(\partial_n s-pn=0\), interpreted variationally; it is not the Neumann
Helmholtz decomposition. Moreover,

\[
 \|F\|_{L^2}+\|p\|_{L^2}+\|s\|_{H^1}+\|g\|_{L^2}
 \le C\|u_h\|_{H(\mathrm{div})}.
 \tag{4}
\]

The pressure is fixed by the all-H1 test space, rather than by an arbitrary
mean-zero convention. Constants in the velocity test space imply the
compatibility condition \(\int F=0\).

## 2. Local regularity at the BDM mesh scale

The required analytic input is the W1,p estimate for the Neumann Stokes
problem on a convex planar domain. Theorem 1.3 of
[Geng and Shen, *Neumann Problems for the Stokes Equations in Convex
Domains*](https://arxiv.org/html/2410.16650v1) supplies it for every finite
1<p<infinity. Its conormal derivative uses the full gradient, exactly as
in (3). We use p=4 on rectangles only.

For each coarse triangle, choose an inner rectangular neighborhood D0
containing it and an outer rectangle D contained in the square, of diameter
comparable to h. At a physical boundary or corner these rectangles meet
that boundary. Choose their sides on coarse grid lines, so D is a union
of complete coarse triangles. A smooth cutoff eta is one on D0 and zero
near the artificial sides of D; its derivatives scale as h^{-j}. Such
neighborhoods can be chosen with uniformly bounded overlap.

The following local estimate holds:

\[
 h^{1/2}\|\nabla s\|_{L^4(D_0)}
 \le C\left(\|\nabla s\|_{L^2(D)}+\|p\|_{L^2(D)}
       +\|g\|_{L^2(D)}+h\|F\|_{L^2(D)}\right).
 \tag{5}
\]

Here are details of the localization, including the inhomogeneous-divergence
point that is not directly part of the cited theorem. Put c=average_D s,
w=eta(s-c), and pi=eta p. Testing (3) with eta times a test vector on D gives

\[
 (\nabla w,\nabla v)-(\pi,\mathrm{div}\,v)
  =(H,v)+(B,\nabla v),
 \quad\mathrm{div}\,w=d,
\]

where

\[
 B=(s-c)\otimes\nabla\eta,\quad
 H=\eta F-\nabla s\,\nabla\eta+p\nabla\eta,\quad
 d=\eta g+\nabla\eta\cdot(s-c).
\]

Scaled Sobolev-Poincare on a rectangle yields
\(\|s-c\|_{L^4(D)}\le C h^{1/2}\|\nabla s\|_{L^2(D)}\).
Since g is constant on each coarse triangle,
\(\|g\|_{L^4(D)}\le C h^{-1/2}\|g\|_{L^2(D)}\).
Thus B and d have L4 norms bounded by h^{-1/2} times the energy/data
terms on the right of (5). Also

\[
 \|H\|_{L^2(D)}\le C\left(\|F\|_{L^2(D)}
       +h^{-1}(\|\nabla s\|_{L^2(D)}+\|p\|_{L^2(D)})\right).
\]

Constant test vectors show that H has zero integral. A componentwise
mean-zero Neumann Poisson solve on D represents (H,v) as (A,grad v), with
\(\|A\|_{L^4(D)}\le C h^{1/2}\|H\|_{L^2(D)}\). This follows by scaling
the usual H2 Neumann estimate and H1 -> L4 embedding on a fixed rectangle.

Take a W1,4 divergence lifting b of d with
\(\|\nabla b\|_{L^4(D)}\le C\|d\|_{L^4(D)}\): a zero-boundary
Bogovskii lifting handles d-average(d), and an affine vector handles the
mean. Now w-b is divergence free and its variational Stokes equation has
an L4 tensor right-hand side B+A-grad b. Applying the cited theorem,
and then adding b back, proves (5). No artificial Dirichlet condition on
the original s, and no assumed H2 regularity of s, is used.

This also proves that s has a continuous representative up to mesh vertices,
including boundary vertices, so its nodal interpolant is well defined.

## 3. Nodal interpolation and the vertex-patch estimate

Let I_h be nodal continuous P1 interpolation, applied componentwise. On a
triangle K, the scaled W1,4 interpolation estimates are

\[
 h^{-1}\|s-I_hs\|_{L^2(K)}+\|\nabla I_hs\|_{L^2(K)}
 \le C h^{1/2}\|\nabla s\|_{L^4(K)}.
\]

Using (5), bounded overlap, and (4) gives

\[
 h^{-1}\|s-I_hs\|_{L^2}+\|I_hs\|_{H^1}
 \le C\|u_h\|_{H(\mathrm{div})}.
 \tag{6}
\]

The canonical BDM1 edge-moment interpolant Pi satisfies, for z in H1(K)^2,

\[
 \|\Pi z\|_{L^2(K)}
 \le C(\|z\|_{L^2(K)}+h\|\nabla z\|_{L^2(K)}).
 \tag{7}
\]

Since I_hs already belongs to BDM1, the local remainder is
r_h=Pi s-I_hs=Pi(s-I_hs). Equations (2), (6), and (7) imply
\(\|r_h\|_{L^2}\le Ch\|u_h\|_{H(\mathrm{div})}\).

On each triangle the vertex split is simply
\(r_{h,v}|_K=\lambda_v(r_h|_K)(v)\). Normal continuity across an edge
is preserved. Finite-dimensional scaling and the inverse estimate give

\[
 \sum_v\|r_{h,v}\|_{H(\mathrm{div})}^2
 \le C\sum_K(1+h_K^{-2})\|r_h\|_{L^2(K)}^2.
 \tag{8}
\]

Combining these estimates proves (1). Constants depend on the fixed domain,
mesh shapes, and analytic estimates, but not on h or on the BDM input.

For these right-isosceles triangles, with h the square side length, an
explicit elementary version of (8) is
\(\sum_v\|r_{h,v}\|_{H(\mathrm{div})}^2
\le(2+24h^{-2})\|r_h\|_{L^2}^2\).
Indeed, for vertex vectors a_i, the local mass norm is
\(|K|(\sum_i|a_i|^2+|\sum_i a_i|^2)/12\), whereas the sum of the
three separate vertex mass norms is \(|K|\sum_i|a_i|^2/6\).
Also \(|\nabla\lambda_i|^2\le2h^{-2}\). This is an upper bound,
not a sharp constant. Controlling div r_h alone would not suffice because
the separate vertex divergences can cancel in their sum.

The potential term is also legitimate, even for an arbitrary admissible
pair with s in H1. Since Rq=u_h-s is H1 on each coarse triangle, q belongs
to H2 on each triangle. Its continuous elementwise representatives agree
along edges and at their endpoints, giving well-defined vertex values.
The vertex-value/edge-average P2 interpolant obeys
\(R\Pi_Qq=u_h-\Pi s\) by edge integration by parts. Scaled H2 interpolation
estimates and the BDM inverse estimate give

\[
 \|\Pi_Qq\|_{H^1}
 \le C\bigl(\|q\|_{H^1}+h|q|_{H^2(\mathcal T_h)}\bigr)
 \le C\bigl(\|q\|_{H^1}+\|u_h\|_{L^2}+h\|s\|_{H^1}\bigr).
 \tag{9}
\]

Here H2(T_h) denotes the broken seminorm. Thus the full original three-term
product norm is uniformly bounded as well, for each fixed positive alpha.
One may additionally choose the discrete potential mean zero to decrease
its norm, but this is not required for stability. The potential transfer
is not being applied to an arbitrary unconstrained H1 scalar.

## 4. A simpler alternative: average only the smooth-space transfer

Let J_h:H1^2 -> continuous P1^2 be an H1-stable local projection with
first-order L2 approximation. A concrete volume-average definition is

\[
 (J_hs)(v)=\frac{1}{|\omega_v|}
       \sum_{K\ni v}\int_K(12\lambda_v-3)s\,dx.
 \tag{10}
\]

On a triangle, (12 lambda_v-3)/|K| is the L2-dual function extracting the
v coefficient of the P1 L2 projection. Formula (10) averages these local
coefficients with area weights and is a projection onto continuous P1.
It applies componentwise, includes boundary vertices, reproduces constants,
and has the local stability/approximation estimates by scaling and
Poincare on vertex neighborhoods.

For r_h=Pi(s-J_hs), (7), the approximation estimate, and H1 stability give

\[
 \sum_v\|r_{h,v}\|_{H(\mathrm{div})}^2\le C\|s\|_{H^1}^2.
 \tag{11}
\]

Unlike the nodal proof, (11) holds for **every** H1 vector s. For the exact
minimizer, (2) proves the desired uniform input bound immediately. For any
admissible approximate minimizer, it bounds the patch energy by its
continuous objective, without needing pointwise convergence.

The same argument for s_a-s_b shows Lipschitz continuity of this patch map
in H1. Thus an actual certified continuous energy-error bound would control
the averaged patch-map error. Successive refinement differences alone are
still not certified bounds against the exact solution.

There is also a genuine convergence conclusion for the existing correction
spaces after this change of transfer. Nested uniform C2 bicubic spline
spaces are dense in H2 on the square. Their Ritz solutions therefore
converge in the correction energy norm. On any fixed BDM mesh the input
space is finite dimensional, so this convergence holds in the corresponding
input-to-continuous-split operator norm. The H1 continuity of the averaged
patch map then proves convergence of its entire input-to-patch map and its
largest generalized eigenvalue to their exact-minimizer values. This is a
convergence theorem, without a certified finite-resolution error or rate.
It does not justify interchanging mesh refinement with an unresolved
continuous solve.

The BDM moment map is already defined on H1 vectors. A matched pair of
bounded commuting projections, such as
[Falk-Winther](https://sites.math.rutgers.edu/~falk/papers/local-cochain-mcomp.pdf),
is another route if one wants transfers defined on the entire de Rham
complex at its natural regularity. It is not necessary merely to obtain
(11). Replacing the potential interpolant independently by an arbitrary
averaging operator would generally lose the commuting identity.

## 5. Why convergence of continuous energy alone was insufficient

The correction equation in the experiment is correct:

\[
 a(\psi,\chi)=(D^2\psi,D^2\chi)+2(\nabla\psi,\nabla\chi)+(\psi,\chi),
\]

with s=s0+R psi and q=q0-psi. Stationarity gives the exact identity

\[
 J(s_\psi,q_\psi)-J(s_*,q_*)=a(\psi-\psi_*,\psi-\psi_*).
\]

Nevertheless, H2 convergence of psi does not control point values of
grad psi in two dimensions. This is a real obstruction, even if the
continuous energy converges to its exact minimum.

For example, fix an interior coarse vertex z and a disk of radius R0
containing no other vertex. There are smooth radial cutoffs f_L equal to
one for r<R0 exp(-L), equal to zero for r>R0, and approximately
log(R0/r)/L between these radii. Smooth the two transitions on comparable
radial scales. For a fixed unit vector e, set

\[
 \delta\psi_L(x)=L^{1/4}\, e\cdot(x-z)\, f_L(|x-z|).
\]

Direct radial integration gives
\(\|\delta\psi_L\|_{H^2}^2=O(L^{-1/2})\), whereas
\(\nabla\delta\psi_L(z)=L^{1/4}e\). Subtracting the scalar mean, if
needed, preserves these conclusions. Fix any nonzero BDM input and perturb
its exact minimizer by
(R delta psi_L,-delta psi_L). The constraint remains exact, and the
objective excess tends to zero. But the nodal P1 vector perturbation grows
without bound at z. Its BDM moment interpolation tends to zero, by (7)
on the fixed coarse mesh. The resulting vertex-patch energy therefore
diverges, while the input u_h is unchanged.

This is **not** a counterexample to (1): those perturbed pairs are not the
exact minimizers. Nor does it prove that the particular spline Galerkin
sequence diverges. It proves that neither a small objective error nor
exact reconstruction alone can certify the nodal patch map. Averaging
removes precisely this discontinuity.

## 6. What the numerical audit says about the slow spline convergence

There is a concrete interior regularity obstruction for the C2 correction
space. Define curl s=partial_x s2-partial_y s1. Compactly supported tests
in the correction equation show that

\[
 (1-\Delta)w=0,\qquad w=\mathrm{curl}\,s-q,
 \qquad (\Delta-1)\psi=\mathrm{curl}\,s_0-q_0-w.
\]

The function w is smooth in the interior. Away from coarse vertices,
transmission across an interior edge e gives

\[
 [D^2\psi]_e=[\mathrm{curl}\,s_0]_e\,n_e\otimes n_e.
\]

Such jumps can occur in all three coarse edge directions. A globally C2
bicubic spline cannot represent them exactly, even when an edge lies on a
spline knot line. Refinement must approximate the jump by a narrow transition.
This is compatible with slow H2 energy convergence; it does not mean the
spline spaces are nondense in H2.

The Opus audit checked the 5x5 maximizing field at relative spline
refinements 8 and 16. At refinement 16, the sums
\(\sum_e[\mathrm{curl}\,s_0]_e^2|e|\) were approximately 20.12, 20.12,
and 21.57 for horizontal, vertical, and diagonal edges. Across diagonal
midpoints, the Hessian difference measured at offsets of two spline cells
agreed with the predicted jumps to about 0.5% in relative Euclidean norm.
At offsets of one quarter cell, the discrepancy was about 43%. These
measurements support a transition spread over roughly a spline cell.

The audit's extrapolations and rate models are heuristic, not certificates
for the limiting patch constants. A triangular C1 macroelement space such
as Hsieh-Clough-Tocher, aligned with the coarse edges, is a concrete future
comparison because it admits Hessian jumps across those edges. It has not
been implemented or tested in this follow-up. Merely repeating tensor knots
on horizontal and vertical grid lines does not resolve the diagonal mismatch.

## 7. Numerical comparison of the two smooth-space transfers

One sequential Opus 5.5 implementation worker assembled (10) with exact
piecewise-polynomial quadrature and compared both transfers using the same
continuous Galerkin solve and the same BDM and patch norm matrices. The
orchestrator did not inspect implementation code. Full numerical results
and checks are in
[the comparison report](agent_artifacts/averaging_comparison.md).

On the 10x10 BDM mesh:

| Spline cells per coarse cell | Nodal maximum | Averaged maximum | Nodal map change | Averaged map change |
|---|---:|---:|---:|---:|
| 4 | 0.691348 | 0.364940 | — | — |
| 8 | 0.454302 | 0.368668 | 0.2814 | 0.01666 |
| 16 | 0.388390 | 0.369386 | 0.1564 | 0.002114 |
| 32 | 0.378613 | 0.369475 | 0.07980 | 0.0002697 |

The maxima are squared patch/input norm ratios. Map changes are unsquared
operator norms from input H(div) to the product of patch H(div) spaces,
relative to the preceding row. At the last refinement the averaged map
change is about 296 times smaller than the nodal change.

On the 5x5 mesh, refinement 32 gives 0.371650 for nodal interpolation and
0.348664 for averaging. The corresponding last map changes are 0.03953 and
0.0001756. Over the tested refinements, averaged map changes decrease by
roughly a factor of eight per doubling, whereas nodal changes decrease by
roughly a factor of two. These are empirical observations, not proven rates.

The maps themselves remain different: at refinement 32 their difference
has operator norm about 0.215 on 10x10 and 0.211 on 5x5. Similar largest
ratios therefore do not mean the decompositions are almost identical for
every input. Averaging gave a smaller maximum in every completed comparison,
but the principal numerical improvement was its much smaller sensitivity
to continuous-solve resolution.

Reconstruction errors were at most 1.4e-10. Continuous P1 reproduction was
accurate to 3.3e-16; independent polynomial and quadrature checks agreed to
roundoff. The nodal columns reproduced the earlier results. These checks
validate the finite-dimensional comparison; they do not certify the limiting
decimal values or demonstrate a convergence rate in the coarse mesh size.
The continuous minimization was identical for the two transfers, so faster
convergence of the averaged patch map does not imply faster convergence of
the underlying correction in H2.

## Scope of the result

The nodal result proved here is for the stated square/grid family and the
stated full-H1 objective with natural boundary conditions. It should not
be silently transferred to a different domain, boundary condition,
continuous decomposition, or independently chosen interpolation pair.
It provides no numerical value for the optimal stability constant.
The previous refinement tables remain approximation diagnostics until a
separate quantitative error bound or convincing independent convergence
study resolves their values.
