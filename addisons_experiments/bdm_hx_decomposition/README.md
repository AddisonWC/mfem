# BDM1 vertex-patch HX decompositions

This experiment uses the unit square divided into a 10×10 grid of squares,
each cut along its lower-left to upper-right diagonal. All vector fields are
represented in the conforming triangular BDM1 space (640 global dimensions).
The script generates 1920×1080 PNGs: four source fields, each with
HX decompositions at α = 1/4, 1, 4; a discrete no-local decomposition; and
an approximate continuous Neumann Helmholtz decomposition.

## Vertex patches

MFEM's BDM1 implementation uses normal values at open edge points. This
experiment uses an equivalent basis whose DOFs are normal values at the two
endpoints of every edge. It is assembled directly from the local six-by-six
interpolation matrix for `[P1]^2`; it does not depend on MFEM's coordinate
ordering. Each endpoint DOF belongs to that endpoint's vertex patch. Its
basis function vanishes at the other two vertices of every incident triangle.
The patches form a direct sum of BDM1 and include boundary vertices. The
script checks this vanishing property numerically.

This follows the endpoint grouping used by the split-vertex
Johnson–Mercier patch construction in `johnson_mercier/common.hpp`. A future
MFEM BDM option could expose this endpoint basis directly; the plot driver
keeps the coordinate change local to the experiment for now.

For BDM coefficient vector `u`, smooth coefficients `s`, potential
coefficients `q`, and patch coefficients `u_v`, the optimizer solves

    minimize  ||s||²_H1 + ||q||²_H1 + α Σ_v ||u_v||²_H(div)
    subject to u = I s + rotgrad q + Σ_v u_v.

Here `s` belongs to continuous vector P1 and `q` to continuous scalar P2.
Both H1 norms include L2 and gradient terms. The local norm includes L2 and
divergence terms. No essential boundary conditions are imposed. `rotgrad q`
means `(∂y q, -∂x q)`. All inner products are integrated exactly for these
polynomial degrees. The Schur complement is the exact inverse additive
preconditioner in the question. The no-local figures solve the same discrete
minimum with all `u_v` omitted; the discrete P1 plus rotgrad P2 map is
surjective for this mesh. This is a discrete comparison, **not** the exact
minimum over infinite-dimensional continuous H1 spaces.

The three source fields are BDM interpolants of two smooth analytic vectors,

    wave:  (sin(πx) cos(πy) + 0.25x, cos(πx) sin(πy) - 0.15y)
    shear: (-0.55 sin(2πy)(0.7+0.3x) + 0.2y,
             0.55 sin(2πx)(0.7+0.3y) + 0.3x),

and a structured field `0.35 wave + rotgrad q_h`. The P2 potential `q_h`
vanishes at every mesh vertex. Its edge midpoint values form a localized
oscillatory pattern with a Gaussian envelope, creating tangential jumps
while preserving normal continuity and giving zero divergence for that part.
The resulting structured field has H(div) norm about 1.12 on this mesh.

The fourth source is a smooth P1 vector (70% wave plus 30% shear at the mesh
vertices) plus `rotgrad q_noise`, where `q_noise` is a continuous P2 potential.
Its values are sampled from four seeded random Fourier bands with integer
frequencies 1–2, 3–4, 5–6, and 7–9. Each band has eight random phases and
directions, and its amplitude is reduced with frequency before P2
interpolation. Bands are then weighted 0.65, 0.48, 0.34, 0.24 in their H1
norms. The combined P2 potential is scaled to have exactly the same assembled
H1 norm as the smooth P1 vector. The random seed is `20260919`.
The [source overview](output/multiscale_noisy_potential_source.png) displays
the scalar potential, the smooth vector, its rotated gradient, and their BDM
sum. Both generating H1 norms are about 2.365.

## Helmholtz figures

The curl-free field is `∇φ`, where `φ` solves the Neumann variational problem

    ∫ ∇φ · ∇ψ = ∫ u · ∇ψ  for every ψ ∈ H1 / constants.

Thus it carries the source's normal boundary flux, and `u - ∇φ` is
divergence-free with zero normal trace in the continuous weak sense. The
figures approximate this continuous projection with scalar P2 elements on a
40×40 grid of squares (3200 triangles). The plotted remainder is the original
BDM field minus that approximate gradient. It is an approximation to the
continuous Helmholtz split, rather than an exactly H(div)-conforming discrete
Helmholtz split.

## Reproduce

From the repository root, with NumPy, SciPy, and Matplotlib installed:

```sh
MPLCONFIGDIR=/tmp/mpl-bdm-hx OPENBLAS_NUM_THREADS=1 \
  python3 addisons_experiments/bdm_hx_decomposition/plot_decompositions.py
MPLCONFIGDIR=/tmp/mpl-bdm-hx OPENBLAS_NUM_THREADS=1 \
  python3 addisons_experiments/bdm_hx_decomposition/plot_noisy_source.py
```

`output/metrics.json` records dimensions, reconstruction errors, and minimum
product norms. The script accepts `--output` and `--helmholtz-resolution`;
the latter must be a multiple of 10 so the fine triangles align with the
BDM mesh for exact load integration.

## Figures

| Source | α = 1/4 | α = 1 | α = 4 | No local | Helmholtz |
|---|---|---|---|---|---|
| Compression wave | [PNG](output/compression_wave_alpha_0.25.png) | [PNG](output/compression_wave_alpha_1.png) | [PNG](output/compression_wave_alpha_4.png) | [PNG](output/compression_wave_no_local.png) | [PNG](output/compression_wave_helmholtz.png) |
| Shear wave | [PNG](output/shear_wave_alpha_0.25.png) | [PNG](output/shear_wave_alpha_1.png) | [PNG](output/shear_wave_alpha_4.png) | [PNG](output/shear_wave_no_local.png) | [PNG](output/shear_wave_helmholtz.png) |
| Localized edge pattern | [PNG](output/localized_edge_pattern_alpha_0.25.png) | [PNG](output/localized_edge_pattern_alpha_1.png) | [PNG](output/localized_edge_pattern_alpha_4.png) | [PNG](output/localized_edge_pattern_no_local.png) | [PNG](output/localized_edge_pattern_helmholtz.png) |
| Multiscale noisy potential | [PNG](output/multiscale_noisy_potential_alpha_0.25.png) | [PNG](output/multiscale_noisy_potential_alpha_1.png) | [PNG](output/multiscale_noisy_potential_alpha_4.png) | [PNG](output/multiscale_noisy_potential_no_local.png) | [PNG](output/multiscale_noisy_potential_helmholtz.png) |

## Constructive HX interpolation experiment

The three additional figures below implement the proof-style construction.
There is a subtlety: the previous Neumann Helmholtz field is a P2 Galerkin
approximation, and directly interpolating its gradient with BDM edge moments
does not give an exact discrete commuting identity. A point-value BDM
interpolant also does not generally commute with divergence.

To retain an **exact continuous identity**, the new driver first obtains the
discrete no-local coefficients `s₀ ∈ [P1]²` and `q₀ ∈ P2`, so
`u_h = s₀ + rotgrad q₀`. It fits a globally smooth Legendre polynomial `ψ`
so that `s = s₀ + rotgrad ψ` approaches the earlier Neumann Helmholtz
curl-free field in L2. It then sets `q = q₀ - ψ`. Thus
`u_h = s + rotgrad q` exactly, with `s ∈ H1²` and `q ∈ H1`. This is a
continuous-level decomposition fitted toward the Helmholtz split, **not the
Helmholtz split itself or its optimal H1 decomposition**. The relative L2
fit error is recorded in the figures and metrics.

`Π_S` is nodal P1 interpolation. `Π_V` uses the two normal-flux P1 moments on
each BDM edge. `Π_Q` uses vertex values and edge averages of the scalar
potential; these define the same P2 space as vertex and midpoint values.
This pair commutes with rotgrad because edge integration by parts uses exactly
the vertex values and edge average. Consequently,

    u_h = (Π_V s - Π_V Π_S s) + Π_V Π_S s + rotgrad Π_Q q.

The first term is partitioned by the endpoint BDM coordinates into vertex
patch functions. Each figure shows the two source components, the two P1
components, the two rotgrad P2 components, and the pointwise sum of local
H(div) densities. It reports squared H1 norms of the P1 and P2 terms,
the sum of squared local H(div) norms, and weighted totals for α = 1/4, 1, 4.
These are norms of this constructive split, not minimized product norms.
At α = 1, the squared norm is 2.39, 19.58, and 2.38 times the corresponding
discrete optimum for the compression, shear, and localized-pattern sources.
The large shear ratio is a property of this chosen fitted continuous split;
the commuting interpolation alone does not make it optimal.

Run with:

```sh
MPLCONFIGDIR=/tmp/mpl-bdm-hx OPENBLAS_NUM_THREADS=1 \
  python3 addisons_experiments/bdm_hx_decomposition/plot_constructive_hx.py
```

| Source | Constructive figure |
|---|---|
| Compression wave | [PNG](output/compression_wave_constructive_hx.png) |
| Shear wave | [PNG](output/shear_wave_constructive_hx.png) |
| Localized edge pattern | [PNG](output/localized_edge_pattern_constructive_hx.png) |
| Multiscale noisy potential | [PNG](output/multiscale_noisy_potential_constructive_hx.png) |

The numeric data are in [constructive_metrics.json](output/constructive_metrics.json).

## H1 × H1 continuous decomposition

The second constructive family addresses the optimal continuous-level
decomposition `u_h = s + rotgrad q`, minimizing
`||s||²_H1 + ||q||²_H1` without local terms. Starting from the exact discrete
no-local pair `(s₀,q₀)`, it writes `s=s₀+rotgrad ψ` and `q=q₀−ψ`. It minimizes
the continuous product norm over globally smooth Legendre polynomials `ψ` of
total degrees 4, 6, and 8. The figure uses degree 8 and reports the decreasing
objective values at all three degrees. This is a finite-dimensional
approximation to the continuous minimum, not a proof of convergence to it.
For any exact split with `s ∈ H1²`, its correction `ψ=q₀−q` has
`rotgrad ψ=s−s₀ ∈ H1²`, so the natural correction space is H2 up to constants.
The constant mode is omitted because the discrete no-local minimizer already
chooses the H1-optimal mean of `q₀`.
The same commuting moment interpolants and endpoint vertex patches give an
exact BDM reconstruction after interpolation.

| Source | H1-product constructive figure |
|---|---|
| Compression wave | [PNG](output/compression_wave_constructive_h1.png) |
| Shear wave | [PNG](output/shear_wave_constructive_h1.png) |
| Localized edge pattern | [PNG](output/localized_edge_pattern_constructive_h1.png) |
| Multiscale noisy potential | [PNG](output/multiscale_noisy_potential_constructive_h1.png) |

The [H1 metrics](output/constructive_h1_metrics.json) include the continuous
and interpolated component norms, the moment-interpolation ratios, and the
degree convergence values. On these four fields, the unsmoothed nodal P1
interpolant and commuting edge moments reconstruct accurately and yield small
local corrections. Point evaluation is not bounded on general H1 functions
in two dimensions, so these experiments do not justify a mesh-independent
bound for nodal interpolation without smoothing.
