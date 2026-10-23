# HX with traction boundary conditions on the stress space

`hx_bc` tests the additive HX preconditioner

\[
 B = R + \Pi_r L_1^{-1}\Pi_r^T + J_r L_2^{-1} J_r^T
\]

for \(a(\sigma,\tau)=(\sigma,\tau)+(\operatorname{div}\sigma,\operatorname{div}\tau)\)
on the stress space \(\Sigma_\Gamma=\{\sigma: \sigma n = 0 \text{ on } \Gamma_s\}\),
where \(\Gamma_s\) is a union of boundary edges. It covers all five pairs on this
branch: JM/HCT, AW/Argyris, HZ/Argyris, HZZZ/Bell, and reduced HZZZ/Bell. The
questions were:

1. Which potential DOFs must be constrained so that the Airy map lands in
   \(\Sigma_\Gamma\)?
2. Is the restricted image the whole div-free subspace of \(\Sigma_\Gamma\)?
3. Is the preconditioner still robust in \(h\)?

**Summary.** With the constraints below, the restricted Airy image is exactly
\(\ker\operatorname{div}\cap\Sigma_\Gamma\). This was checked by dense ranks on
every configuration up to 8×8. HX stays \(h\)-robust for every \(\Gamma_s\)
tried. The condition number grows from 5.1 (no BC) to 9–20 depending on
\(\Gamma_s\), or 22–31 PCG iterations instead of 15. The constant does not
depend on the element pair. Two choices matter:

* **Argyris/Bell vertex Hessians:** impose \(H t = 0\) for each incident
  \(\Gamma_s\) tangent; do not zero the whole Hessian.
* **H1 auxiliary space:** impose only \(\sigma n = 0\), not \(\sigma = 0\).

## DOF analysis

Write \(R\) for rotation by +90° and \(n = R t\) on an edge with unit tangent
\(t\). Then

\[
 \operatorname{Airy}(\phi) = \begin{pmatrix}\phi_{yy}&-\phi_{xy}\\-\phi_{xy}&\phi_{xx}\end{pmatrix},
 \qquad \operatorname{Airy}(\phi)\,n = R\,\partial_t \nabla\phi .
\]

So \(\sigma n=0\) on \(\Gamma_s\) means \(\nabla\phi\) is constant along each
connected component \(\Gamma_s^j\), and \(\phi\) agrees there with an affine
function \(a_j + b_j\cdot x\). The affine kernel of Airy removes
\((a_1,b_1)\): \(\phi\) is **clamped** on \(\Gamma_s^1\). On each further
component, \(\phi\) has the trace of a *free* affine function. That gives 3
global parameters per extra component (the cohomology).

### Stress spaces

| DOF | on \(\Gamma_s\) |
|---|---|
| edge traction moments (AW/HZ: 4, HZZZ: 3, JM: 4) | all zero |
| AW/HZ/HZZZ vertex (xx,xy,yy), one incident direction \(n\) | \(\sigma n=0\): only \(\sigma = s\,t\otimes t\) remains (1 DOF) |
| same, at a corner of \(\Gamma_s\) (two directions) | \(\sigma(v)=0\) |
| JM | no vertex DOFs |

The vertex DOFs are Cartesian, so on an oblique boundary this is a change of
basis, not a DOF elimination. On an axis-aligned side \(y=c\), eliminate
\((xy,yy)\); on \(x=c\), eliminate \((xx,xy)\). At a point where
\(\Gamma_s\) ends, only the \(\Gamma_s\) edge's normal counts.

### Potential spaces (the Argyris question)

| DOF (Argyris: v, ∂x, ∂y, xx, xy, yy per vertex + edge ∂n) | on \(\Gamma_s^1\) |
|---|---|
| value, gradient | zero |
| Hessian, one incident tangent \(t\) | \(H t = 0\): \(\phi_{tt}=\phi_{nt}=0\), \(\phi_{nn}\) **free** |
| Hessian at a corner of \(\Gamma_s\) | \(H = 0\) |
| edge midpoint normal derivative | zero |

These are exactly the clamped conditions on each edge. \(\phi|_e\) is quintic
and fixed by \((\phi,\phi_t,\phi_{tt})\) at the endpoints. \(\partial_n\phi|_e\)
is quartic and fixed by \((\partial_n\phi,\partial_t\partial_n\phi)\) at the
endpoints plus the midpoint DOF. On a side \(y=c\), eliminate \((xx,xy)\); on
\(x=c\), eliminate \((xy,yy)\). Because \(\operatorname{cof}\) maps \(n\otimes n\)
to \(t\otimes t\), the free \(\phi_{nn}\) is exactly the free vertex
\(\sigma_{tt}\). The Hessian constraint is the vertex stress constraint pulled
back.

* **Bell:** the same vertex conditions, and it has no edge DOFs. Its cubic
  \(\partial_n\phi|_e\) is fixed by the endpoint data.
* **HCT:** value and gradient zero at the vertices, plus the edge DOF. There
  are no Hessian DOFs, so there is nothing to decide.
* **Components \(j\ge2\):** the same constrained DOFs take the values of
  \(1\), \(x\), \(y\), as three extra columns. The vertex jet is exact (value,
  gradient, zero Hessian); edge DOFs come from a per-element fit.
* **Empty \(\Gamma_s\):** the existing point gauge (value and gradient at
  vertex 0).

### H1 auxiliary space

Continuous P1 tensors with \(\sigma n=0\) at \(\Gamma_s\) vertices (`-h1bc
normal`), or \(\sigma=0\) there (`-h1bc full`). Either way the canonical
interpolant lands in \(\Sigma_\Gamma\): edge tractions are linear and vanish
at both endpoints.

## Implementation

Each space gets a sparse \(P\) whose columns span the constrained subspace.
These are identity columns, orthonormal vertex null-space columns, and the
cohomology columns. Then \(A_r=P_s^TAP_s\),
\(\Pi_r=P_s^T\Pi P_{h1}\), \(J_r=P_s^T\,\mathrm{Airy}\,P_\phi\), and the
auxiliary operators are \(P^TLP\) (weighted vector H1, Hessian form). \(P_s\)
has orthonormal columns. `leak` = \(\max|(I-P_sP_s^T)M|/\max|M|\) checks that
\(M=\Pi P_{h1}\) or \(\mathrm{Airy}\,P_\phi\) lands in \(\Sigma_\Gamma\).

Other details:

* **Exactness:** for restricted sizes up to 2500, `ker_div` is the number of
  generalized eigenvalues of (div-div, mass) below 1e-8. `rank_airy` is the
  numerical rank of \(J_r\).
* **Smoother:** the ex43 vertex patches, mapped to restricted coordinates and
  scaled by 0.33.
* **Auxiliary solves:** exact (UMFPACK).
* **PCG:** MFEM PCG to relative tolerance 1e-8 in the preconditioned norm,
  with a seeded random load and a zero initial guess.
* **κ:** Ritz extremes of \(BA\) by A-inner-product Lanczos with full
  reorthogonalization. It stops when both extremes change by less than 1e-6
  relative over 5 steps, or at 300 steps (200 at level 6). Ritz values are
  inner bounds, so κ is a converging lower estimate.

Mesh: the unit square, \(N\times N\) squares split into two triangles,
\(N = 2^{\ell+1}\). \(\Gamma_s\) codes:

| code | \(\Gamma_s\) |
|---|---|
| `b` | bottom |
| `h` | lower half of left (ends mid-side) |
| `lr` | left + right (2 components) |
| `bl` | bottom + left (one corner) |
| `blt` | U shape |
| `a` | whole boundary |

## Results

Entries are `PCG iterations / κ`. Full data: [`results/hx_bc.csv`](../results/hx_bc.csv).

### Correct constraints, patch smoother, AW

| level | h | none | h | b | lr | bl | a | blt |
|---|---|---|---|---|---|---|---|---|
| 0 | 1/2 | 12 / 2.7 | 13 / 3.0 | 13 / 3.5 | 15 / 4.5 | 14 / 4.2 | 15 / 5.1 | 17 / 11.2 |
| 1 | 1/4 | 13 / 3.5 | 15 / 4.0 | 17 / 4.9 | 18 / 8.8 | 18 / 6.9 | 22 / 8.5 | 23 / 17.2 |
| 2 | 1/8 | 15 / 4.5 | 18 / 5.7 | 19 / 6.9 | 21 / 10.9 | 21 / 9.3 | 24 / 10.8 | 27 / 19.5 |
| 3 | 1/16 | 15 / 5.0 | 20 / 7.6 | 21 / 8.3 | 23 / 11.4 | 23 / 10.9 | 25 / 11.9 | 28 / 20.1 |
| 4 | 1/32 | 15 / 5.1 | 22 / 8.6 | 22 / 8.8 | 23 / 11.5 | 24 / 11.9 | 25 / 12.5 | 29 / 20.3 |
| 5 | 1/64 | 15 / 5.1 | 22 / 8.9 | 22 / 9.0 | 23 / 11.6 | 23 / 12.6 | 24 / 12.8 | 29 / 20.3 |
| 6 | 1/128 | 14 / 5.1 |  |  | 23 / 11.6 | 24 / 13.1 | 24 / 13.0 | 29 / 20.4 |

The growth at the coarse levels is pre-asymptotic: boundary vertices are a
large fraction of the mesh there. `bl` is the slowest to settle (level-to-level
increments 1.0, 0.68, 0.47), consistent with a limit near 14. The constant
depends on the shape of \(\Gamma_s\) and is not monotone in its length: the U
shape is worse than the whole boundary. I have not found the mode responsible.

### All element pairs, level 5

| Γ_s | jm | aw | hz | hzzz | hzzzr |
|---|---|---|---|---|---|
| none | 16 / 5.2 | 15 / 5.1 | 15 / 5.1 | 15 / 5.1 | 15 / 5.1 |
| h | 23 / 8.9 | 22 / 8.9 | 22 / 8.9 | 22 / 8.9 | 22 / 8.9 |
| b | 22 / 8.9 | 22 / 9.0 | 22 / 9.0 | 22 / 9.0 | 22 / 9.0 |
| lr | 25 / 11.6 | 23 / 11.6 | 23 / 11.6 | 23 / 11.6 | 23 / 11.6 |
| bl | 27 / 12.5 | 23 / 12.6 | 23 / 12.6 | 23 / 12.6 | 23 / 12.6 |
| a | 27 / 12.6 | 24 / 12.8 | 24 / 12.8 | 24 / 12.8 | 24 / 12.8 |
| blt | 31 / 20.3 | 29 / 20.3 | 29 / 20.3 | 29 / 20.3 | 29 / 20.3 |

### Exactness (levels 0, 1, 2)

| element | Γ_s | variant | comps | dim ker div | rank \(J_r\) | max leak |
|---|---|---|---|---|---|---|
| jm | bl | correct | 1 | 24, 96, 384 | 24, 96, 384 | 7e-15 |
| aw | bl | correct | 1 | 40, 152, 592 | 40, 152, 592 | 5e-12 |
| hzzz | bl | correct | 1 | 28, 104, 400 | 28, 104, 400 | 2e-12 |
| aw | a | correct | 1 | 18, 106, 498 | 18, 106, 498 | 5e-12 |
| aw | lr | correct | 2 | 39, 151, 591 | 39, 151, 591 | 3e-12 |
| aw | lr | no cohomology | 2 | 39, 151, 591 | 36, 148, 588 | 3e-12 |
| aw | bl | Hessian zero | 1 | 40, 152, 592 | 36, 144, 576 | 3e-12 |
| aw | bl | Hessian free | 1 | 40, 152, 592 | 44, 160, 606 | 0.25 |
| hzzz | bl | Hessian free | 1 | 28, 104, 400 | 36, 120, 432 | 0.25 |

Every correct configuration in the sweep is exact: 5 elements × 7 \(\Gamma_s\),
rotated and jittered meshes included.

### Variants

| element | Γ_s | variant | L0 | L1 | L2 | L3 | L4 | L5 | L6 |
|---|---|---|---|---|---|---|---|---|---|
| aw | bl | correct | 14 / 4.2 | 18 / 6.9 | 21 / 9.3 | 23 / 10.9 | 24 / 11.9 | 23 / 12.6 | 24 / 13.1 |
| aw | bl | Hessian zero | 21 / 9.8 | 25 / 10.8 | 27 / 11.4 | 28 / 11.7 | 28 / 11.9 | 27 / 12.6 |  |
| aw | bl | Hessian free | 18 / 10.2 | 28 / 33.1 | 50 / 142 | 98 / 619 | 204 / 2666 | 402 / 11228 |  |
| hzzz | bl | Hessian zero | 22 / 13.5 | 30 / 15.6 | 33 / 17.0 | 34 / 17.6 | 34 / 17.9 | 33 / 18.0 |  |
| hzzz | bl | Hessian free | 27 / 86 | 57 / 541 | 102 / 2866 | 167 / 13280 | 230 / 57985 | 269 / —\* |  |
| aw | bl | H1 σ=0 | 16 / 8.0 | 22 / 17.8 | 28 / 29.5 | 34 / 40.0 | 37 / 48.8 | 39 / 55.9 | 40 / 61.5 |
| jm | bl | H1 σ=0 | 16 / 8.3 | 23 / 18.0 | 29 / 29.6 | 35 / 40.0 | 40 / 48.7 | 42 / 55.8 |  |
| aw | a | H1 σ=0 | 19 / 6.7 | 28 / 17.3 | 37 / 29.9 | 43 / 41.1 | 45 / 49.6 | 42 / 56.1 |  |
| aw | lr | no cohomology | 19 / 15.4 | 23 / 25.1 | 28 / 28.6 | 29 / 29.6 | 30 / 29.9 | 30 / 30.0 | 29 / 30.1 |
| aw | a | rotated 30° | 15 / 5.1 | 21 / 8.5 | 25 / 10.8 | 24 / 11.9 | 25 / 12.5 | 25 / 12.8 |  |
| hzzz | bl | rotated 30° | 14 / 4.2 | 18 / 6.9 | 22 / 9.3 | 24 / 10.9 | 25 / 11.9 | 26 / 12.6 |  |
| jm | a | jitter 0.2h | 17 / 4.9 | 23 / 8.6 | 26 / 10.7 | 27 / 11.6 | 27 / 12.2 | 27 / 12.6 |  |
| aw | none | Jacobi | 24 / 13.1 | 38 / 33.9 | 51 / 64.6 | 60 / 96.8 | 58 / 111 | 56 / 115 |  |
| aw | a | Jacobi | 46 / 127 | 66 / 105 | 70 / 112 | 72 / 115 | 68 / 116 | 65 / 116 |  |
| aw | none | sym. GS | 16 / 4.3 | 21 / 8.5 | 26 / 15.2 | 28 / 20.6 | 27 / 22.8 | 26 / 23.5 |  |
| aw | a | sym. GS | 30 / 31.1 | 36 / 24.1 | 35 / 23.3 | 34 / 23.4 | 32 / 23.5 | 30 / 23.5 |  |

\* Lanczos did not converge in 300 steps; the CSV value (4.1e4) is meaningless.

Observations:

* **Free Hessian** (clamp only value, gradient, and edge DOF, as for HCT): the
  Airy image leaves \(\Sigma_\Gamma\). After projection, \(P_s^TJ\) is no longer
  div-free, and \(\lambda_{\max}\) grows like \(h^{-2}\). This is the one
  catastrophic mistake. The `hx_bc_free_hessian_fails` test pins it down.
* **Zero Hessian** (\(\phi\) "fully clamped" at \(\Gamma_s\) vertices): the
  image stays in \(\Sigma_\Gamma\) but misses one \(\sigma_{tt}\) mode per
  non-corner \(\Gamma_s\) vertex. These are local, the patch smoother catches
  them, and the result is still robust. It costs 3–4 iterations for AW and
  about 10 for HZZZ (κ 18 vs 12.6). Keep \(\phi_{nn}\) free.
* **H1 with σ=0 at Γ_s vertices:** κ keeps growing (30 → 40 → 49 → 56 → 61,
  increments 10.5, 8.8, 7.1, 5.7). At best it saturates far higher; it may grow
  logarithmically. My guess is that the boundary-layer \(\sigma_{tt}\) with
  nonzero divergence then has no stable auxiliary representation. This is
  unverified. With \(\sigma n=0\) only, the same configurations are bounded
  near 13.
* **No cohomology columns:** the 3-dim gap costs a bounded factor (κ 11.6 →
  30, 23 → 30 iterations) but no \(h\)-growth. The missing div-free fields are
  smooth, so the H1 space represents them stably. If the global affine columns
  are inconvenient (multigrid, parallel), dropping them is a reasonable
  trade-off.
* **Rotation:** the rotated square reproduces the axis-aligned numbers to all
  printed digits. This confirms the Cartesian-vertex change of basis. Interior
  jitter changes the constants only slightly.
* **Smoothers:** Jacobi and symmetric GS converge to the same asymptotic κ with
  and without BC (≈116 and ≈23.5 for AW). The BC penalty is specific to how
  well the patch smoother does on the interior problem.

## Limitations and next steps

* **Serial, direct solves:** the runs are serial and assembled, with exact
  auxiliary solves. The `-gmg` path of `ex43_hx` would need restricted
  hierarchies. Argyris and P1 are nested and clamped coarse functions stay
  clamped, so \(P_f^T T P_c\) should work. HCT and Bell use DOF-interpolation
  transfer, which needs checking. Cohomology columns would then need special
  handling or omission.
* **Domain:** only the unit square and straight sides. Curved boundaries are
  outside the affine-element assumptions anyway.
* **Mixed BCs:** displacement (natural) BCs on the remaining boundary are
  untested beyond leaving it free. A re-entrant corner of \(\Gamma_s\) is also
  untested.

## Reproduction

Build the project as in the [README](../../README.md). `hx_bc` is built when
CMake finds LAPACK. From `addisons_experiments/build`:

```sh
sh ../johnson_mercier/hx_bc_sweep.sh 6 > hx_bc.csv      # about 40 minutes
OPENBLAS_NUM_THREADS=1 ./hx_bc -e aw -bc bl -rmax 5     # one configuration
./hx_bc -e hzzz -bc h -rmax 1 -check                    # exit 3 on any failure
```

Options:

| option | meaning |
|---|---|
| `-e` | element: `jm`, `aw`, `hz`, `hzzz`, `hzzzr` |
| `-bc` | traction sides: any of `b`,`r`,`t`,`l`,`h`, or `a`, `none` |
| `-h1bc` | `normal` or `full` |
| `-hess` | `tangent` (correct), `free`, `zero` |
| `-no-cohom` | drop the cohomology columns |
| `-s` | `patch`, `jacobi`, `gs` |
| `-rot` | rotation in degrees |
| `-jitter` | interior vertex jitter, as a fraction of h |
| `-rmin`/`-rmax` | refinement range |
| `-lanczos` | maximum Lanczos steps |
| `-exact-max` | largest restricted size for the dense checks |

Provenance:

* **Source:** commit `8e633997ab` plus the then-uncommitted `hx_bc.cpp`,
  `hx_bc_sweep.sh`, and CMake changes. The level 0–5 rows ran before two final
  edits: adding `-check`, and exact vertex jets in the cohomology columns. The
  only rows those affect (`lr` with cohomology) were regenerated with the final
  code; the other rows are unaffected.
* **Build:** MFEM Release build (`build-hx`), GCC 16.2.1, SuiteSparse UMFPACK,
  OpenBLAS 0.3.29 with `OPENBLAS_NUM_THREADS=1`.
* **Hardware:** AMD Ryzen 7 5825U, 6 concurrent jobs; level 6 ran serially
  (about 3.5 GB each).
