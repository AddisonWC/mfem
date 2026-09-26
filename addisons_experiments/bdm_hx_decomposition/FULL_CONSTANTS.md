# Full decomposition constants, alpha = 1

All constants below use squared norms and the full H(div) input norm.
They are not the patch-only constants in the earlier interpolation table.

## Definitions

Let A(u)=||u||^2_H(div). Define

\[
 C_{\rm opt,HX}(h)=\sup_{0\ne u_h\in V_h}
 \frac{\min_{u_h=s_h+R q_h+\sum_v r_v}
       (\|s_h\|_{H^1}^2+\|q_h\|_{H^1}^2+
                   \sum_v\|r_v\|_{H(\mathrm{div})}^2)}{A(u_h)}.
\]

Here s_h is continuous vector P1, q_h is continuous scalar P2, and the
r_v belong to the specified endpoint vertex spaces, with no essential
boundary conditions. This is the optimal discrete HX constant.

For the exact continuous minimizer u_h=s+Rq, the explicit averaged
construction instead uses

\[
 s_h=J_hs,\quad q_h=\Pi_Qq,\quad
 \sum_vr_v=\Pi_Vs-J_hs.
\]

Its full constant C_avg(h) is the supremum of that sum of three energies
divided by A(u_h). The energy sum must be maximized as one quadratic form;
adding the separate component maxima generally overestimates it. A scalar
constant may optionally be subtracted from q_h to reduce its H1 norm
without changing reconstruction; numerical reports must specify whether
this normalization is made.

For comparison define C_cont(h) by replacing the discrete minimization
above by min(||s||^2_H1+||q||^2_H1) over the continuous constraint u_h=s+Rq.
The finite spline approximations to this third constant are upper
approximations and may substantially exceed the exact value.

## Exact continuous benchmark: C_cont(h) tends to 1

In fact the optimal continuous constant over all of H(div) on the square
is exactly 1. The proof uses the full H1 norm, including the potential mass.

Given u in H(div), let g=div u and solve

\[
 \Delta\phi=g,\qquad \phi|_{\partial\Omega}=0.
\]

Set a=grad phi. On the square, phi is H2 and
\(\|D^2\phi\|_{L^2}^2=\|g\|_{L^2}^2\). The vector w=u-a is
divergence free and L2-orthogonal to a. Write w=Rq0 with q0 in H1 and
mean zero. Solve the scalar Neumann resolvent

\[
 p-\Delta p=q_0,\qquad \partial_n p=0.
\]

Take s=a+Rp and q=q0-p. This pair reconstructs u. Both a and Rp have
zero tangential trace. On the rectangle the full-gradient div-curl
identity holds for H1 vectors with zero tangential trace. It implies

\[
 (a,Rp)_{H^1}=0,\qquad
 \|a\|_{H^1}^2=\|a\|_{L^2}^2+\|g\|_{L^2}^2,
 \qquad\|Rp\|_{H^1}^2=\|\nabla p\|_{L^2}^2+\|\Delta p\|_{L^2}^2.
\]

These identities can also be checked directly with the sine basis for
the Dirichlet problem and the cosine basis for the Neumann problem.
The resolvent gives q=-Delta p in H1 and

\[
 \|Rp\|_{H^1}^2+\|q_0-p\|_{H^1}^2
 =\|\nabla p\|^2+2\|\Delta p\|^2+\|\nabla\Delta p\|^2
 =\|\nabla q_0\|^2=\|w\|^2.
\]

Consequently this admissible pair has energy **exactly**
\(\|u\|_{H(\mathrm{div})}^2\). It need not itself minimize that energy,
but it proves the optimal continuous energy is at most the input energy.

To attain equality, choose any nonzero phi in C_c^infinity(Omega) and
u=grad phi. The pair (s,q)=(u,0) is stationary for every admissible
variation (R chi,-chi), chi in H2, since
\((\nabla\phi,R\chi)_{H^1}=0\) by integration by parts with compact
support. Strict convexity makes it the minimizer. Its H1 energy equals
the H(div) energy of u. Hence the continuous constant is exactly 1.

Finally approximate this fixed smooth u by its BDM interpolants u_h in
H(div). The continuous minimal-preimage map is linear and bounded by the
energy estimate just proved, so its energy converges on these inputs.
It follows that

\[
 C_{\rm cont}(h)\le1\quad\hbox{for every }h,
 \qquad\lim_{h\to0}C_{\rm cont}(h)=1.
\]

Thus finite-spline values above 1 are not estimates of a continuous
constant greater than 1: they contain continuous discretization error.

## A lower bound for the discrete HX limit

The continuous result does not by itself identify C_opt,HX. Nevertheless,

\[
 \liminf_{h\to0}C_{\rm opt,HX}(h)\ge1,
 \qquad C_{\rm avg}(h)\ge C_{\rm opt,HX}(h).
\]

For the first statement use u=grad phi with phi smooth and compactly
supported, as above, and u_h converging to u in H(div). Test any discrete
decomposition of u_h against u in the H(div) inner product. The potential
term vanishes. The smooth term equals (s_h,u)_H1. An interior vertex term
r_v has zero total divergence on its support, so

\[
 (r_v,u)_{H(\mathrm{div})}
 =(\mathrm{div}\,r_v,\Delta\phi-\phi)
 =(\mathrm{div}\,r_v,
       (\Delta\phi-\phi)-(\Delta\phi-\phi)(v)).
\]

Boundary patches do not meet the support of phi once h is small. Scaling
and bounded overlap therefore bound the sum of these local terms by
\(h K_\phi (\sum_v\|r_v\|_{H(\mathrm{div})}^2)^{1/2}\), for a fixed
K_phi independent of h. Cauchy-Schwarz across the product space gives

\[
 E_{\rm opt,HX}(u_h)
 \ge\frac{(u_h,u)_{H(\mathrm{div})}^2}
 {\|u\|_{H^1}^2+h^2K_\phi^2}.
\]

Divide by A(u_h) and pass to the limit to obtain the lower bound 1.
No matching upper bound of 1 for the discrete HX optimum is asserted here.
The inequality C_avg >= C_opt,HX holds because the explicit construction
is one admissible discrete decomposition.

## Completed numerical mesh study (September 25, 2026)

The sequential Opus 5.5 worker completed n=30 at spline refinements 8 and 16.
The selected [final numerical report](results/full_constants/full_constants_final.md)
and [compact data and fits](results/full_constants/full_constants_summary.json)
include reproduction commands and provenance. These files are outside ignored
routine output. They have not been committed.

| n | Optimal discrete HX | Averaged, r=8 | Averaged, r=16 | Averaged, r=32 |
|---|---:|---:|---:|---:|
| 5 | 0.997165842 | 1.070342557 | 1.069427462 | 1.069311349 |
| 10 | 0.999378012 | 1.144159892 | 1.143055714 | 1.142915510 |
| 15 | 0.999737259 | 1.159968698 | 1.158828209 | 1.158683368 |
| 20 | 0.999855796 | 1.165510575 | 1.164355940 | — |
| 30 | 0.999937301 | 1.169434689 | 1.168269179 | — |

All values use full squared norms and alpha=1. The averaged values use the
finite spline minimizer; the definition of C_avg above uses the exact continuous
minimizer. Potential mean subtraction changes the reported maxima only at
rounding level.

### Interpretation of the limiting estimates

The numerical estimate for the optimal discrete HX limit is **1**. Larger-mesh
fits give intercepts within about 2e-5 of 1. The proof above establishes only
liminf >= 1 for this discrete constant; no matching upper bound is claimed.

For the explicit averaged construction, the data suggest **about 1.17**
(approximately **1.171** if an additional digit is useful). Fits of the r=16
series at n=10,15,20,30 using L+a/n^2, L+a/n^2+b/n^3 and L+a/n+b/n^2 give
intercepts 1.171436, 1.171412 and 1.171367. Restricting to n=15,20,30 gives
1.171425, 1.171366 and 1.171252; the three-parameter fits on three points
interpolate exactly. Including n=5 widens the fit spread to 1.1701–1.1739.
These variations are diagnostics, not confidence intervals or rigorous bounds.

Spline error remains distinct from mesh extrapolation. At n=30, refinement
8→16 lowers the averaged maximum by 0.001166 and changes the full product map
by 0.02091 in operator norm. The scalar change resembles that on the smaller
meshes, but the map change grows with n. A speculative continuation of the
refinement trend seen at n<=15 would lower the r=16 fit intercept by about
0.00017; that correction has not been verified on the larger meshes.

In particular, this study does not justify exchanging h→0 and spline-resolution
limits. The spline continuous objective at n=30 is still 41.3 at r=16, compared
with the exact bound C_cont(h)<=1. Averaging makes the output constant much less
sensitive to that error, but neither successive differences nor small algebraic
residuals certify its error against the exact continuous-minimizer construction.
The exact continuous benchmark, the empirical discrete limit, and the empirical
averaged limit must remain separate claims.

### Validation and stopping point

At n=30 the optimal discrete reconstruction error is 1.6e-15 relative, and the
averaged reconstruction error is below 6e-9. Generalized-eigen residuals are
at most 3.4e-12. The reported energy identities and numerical form dominance
C_opt,HX <= C_avg hold; details and rounding tolerances are in the report.
The parent inspected the numerical report and mathematical definitions, not
implementation code, and did not rerun the worker's checks.

The n=30 computation completed successfully; its logged run took about
48 minutes and peaked at 11.4 GB of memory. No further sweeps are planned for
this wrap-up. Obtaining a certified exact-limit value would be a separate
research task, rather than another decimal-place extrapolation.
