# Johnson–Mercier multigrid with split patches and HCT Jacobi

The split smoother now adds Jacobi corrections along the Airy images of the
HCT potential basis to the existing split vertex patches. This substantially
reduces iteration growth, but does **not** demonstrate refinement-independent
convergence. On the original refinement range, macro patches remain faster.

`mg_compare` solves the same serial, assembled stress problem as `ex43` and
`hx_compare`, with no essential stress boundary conditions:

\[
 A(\sigma,\tau)=(\sigma,\tau)+(\mathrm{div}\,\sigma,\mathrm{div}\,\tau).
\]

The meshes are connected, straight, conforming triangulations. Every hierarchy
prefix uses a symmetric V-cycle, equal pre/post smoothing counts, and a shared
fixed direct coarse solve.

## Smoother and unchanged transfer

The hierarchy is constructed once in the original moment JM collection. All
variants borrow identical prolongation objects and rediscretized level
operators; restriction is the transpose of prolongation. JM macro refinements
are generally nonnested: this is moment interpolation, not exact inclusion or
Galerkin coarsening.

Let `B` map split coefficients to moments, and let `Av = B^T A B`. The existing
split patch inverse in moment coordinates is

\[
 R_p=B\left[\sum_p I_p(I_p^T A_v I_p)^{-1}I_p^T\right]B^T.
\]

Let `C` be the HCT-to-moment-JM Airy interpolation matrix. For a scalar potential
\(\phi\), its Airy stress is
\((\phi_{yy},-\phi_{xy};-\phi_{xy},\phi_{xx})\), which has zero divergence.
The new `split` smoother is the additive two-term operator

\[
 S r=\omega_p R_p r+\omega_a C D^{-1}C^T r,
 \qquad D_{ii}=(C e_i)^T A(C e_i).
\]

Both terms act on the same residual. Thus the second term is a sum of
one-dimensional corrections along the images of **all** HCT basis functions.
It is a Hiptmair-style potential-space Jacobi smoother, not an exact auxiliary
solve. Since `div C = 0`, `D` is assembled as the diagonal of the HCT Hessian
form, avoiding cancellation from the div-div part of `A`. The affine potential
kernel requires no gauge: individual basis images have positive energy, and
only the diagonal is inverted. The HCT space lives on each level's macro mesh;
its piecewise cubics use the same Alfeld subdivision as JM.

| Name | Smoothing correction |
|---|---|
| `macro` | Original vertex stars in moment coordinates |
| `macro-vertex` | Same stars factored in split coordinates; basis-invariance control |
| `split-patches` | Previous split smoother: edge-endpoint traction patches and separate three-DOF barycenter blocks |
| `split` | Split patches plus HCT Airy-image Jacobi |

`-damping` sets \(\omega_p\), default 0.33; `-airy-damping` sets
\(\omega_a\), default 0.05. All patch weights must lie in `(0,0.5)`.
For `split` (including `all`), the driver additionally requires positive Airy
weight and

\[
 4\omega_p+12\omega_a<2.
\]

At most four split patch subspaces and twelve HCT basis images contribute on
one macro element. Elementwise Cauchy–Schwarz therefore bounds the largest
eigenvalue of the weighted additive correction times `A` by this sum. The
default is 1.92, below 2, giving a conservative sufficient condition for positive
symmetric smoothing in PCG. These weights were not optimized. The previous
`-damping 0.49` setting needs a smaller Airy weight if used with `split`.

## Reproduction and provenance

Measurements below were made on 2026-09-15 from base revision
`d4f4a67a69658c1ffcb3610700655c41e47b5e72` plus the accompanying uncommitted
changes to the experiment's shared solver, multigrid driver, and tests.
Both MFEM (`build-hx`) and the standalone experiment use Release builds;
compiler: GCC 16.2.1 20260819 (Red Hat 16.2.1-2); direct backend: SuiteSparse
UMFPACK; CPU: AMD Ryzen 7 5825U (8 cores, 16 logical CPUs).
`OPENBLAS_NUM_THREADS=1` was used throughout. Timings are single-run wall times
on a shared machine, not controlled performance estimates.

From the repository root, with MFEM already built as described in the
[README](../../README.md):

```sh
cmake -S addisons_experiments -B addisons_experiments/build \
  -DMFEM_DIR="$PWD/build-hx" -DCMAKE_BUILD_TYPE=Release
cmake --build addisons_experiments/build -j 4
ctest --test-dir addisons_experiments/build --output-on-failure
mkdir -p addisons_experiments/output
OPENBLAS_NUM_THREADS=1 addisons_experiments/build/mg_compare -r 5 \
  > addisons_experiments/output/mg_triangle.csv
OPENBLAS_NUM_THREADS=1 addisons_experiments/build/mg_compare \
  -m addisons_experiments/build/data/inline-tri.mesh -r 3 \
  > addisons_experiments/output/mg_square.csv
OPENBLAS_NUM_THREADS=1 addisons_experiments/build/mg_compare -r 5 -s 2 \
  > addisons_experiments/output/mg_two_sweeps.csv
OPENBLAS_NUM_THREADS=1 addisons_experiments/build/mg_compare -r 6 -smoother split \
  > addisons_experiments/output/mg_triangle_deeper.csv
```

All four benchmark commands above return status 3: the first three because
the patch-only baseline fails tolerance, and the fourth because the updated
split smoother misses `1e-8` on the extra refinement. Use `-smoother split` to run only
the updated preconditioner. `-cr` refines before constructing the coarse grid;
`-r` adds hierarchy levels; `-s` sets steps on each side of the coarse correction.
Without SuiteSparse, the coarse solver uses dense LU and limits its matrix to
2000 unknowns. HCT Jacobi uses no auxiliary direct solve.

Every solve starts from zero with the same seeded random moment load at a given
depth. Except for the explicitly relaxed run below, PCG must meet a recomputed
unpreconditioned relative residual of `1e-8`,
with a limit of 2000 iterations. CSV `setup_s` includes patch factorization,
HCT map/diagonal assembly for `split`, and V-cycle construction. It excludes
shared meshes, transfers, level operators, basis maps, and coarse factorization.
`solve_s` measures one solve. The new `airy_damping` column is zero for variants
without the HCT term.

## Results

[mg_hct_results.csv](../results/mg_hct_results.csv) records the current runs,
with added `dataset` and `tolerance` columns identifying the command and target.
There are 85 rows, including retained failures. The older
[mg_results.csv](../results/mg_results.csv) remains a historical snapshot;
its `split` label means the present `split-patches` baseline.

Reference triangle, one pre/post step:

| JM unknowns | Macro (either basis) | Split patches only | Split + HCT Jacobi |
|---:|---:|---:|---:|
| 15 | 1 | 1 | 1 |
| 48 | 13 | 41 | 37 |
| 168 | 20 | 141 | 83 |
| 624 | 27 | 362 | 110 |
| 2,400 | 32 | 856 | 134 |
| 9,408 | 34 | 2000 (failed) | 152 |

At 9,408 unknowns, the new smoother reaches relative residual `9.49e-9`;
the baseline reaches only `2.12e-8`. Solve times are approximately 0.27 s for
the new smoother, 2.50 s for the failed baseline, and 0.083 s for macro patches.
Setup is approximately 0.072 s, 0.0023 s, and 0.010 s, respectively: assembling
the HCT correction increases setup cost.

Square (`inline-tri.mesh`), one pre/post step:

| JM unknowns | Macro (either basis) | Split patches only | Split + HCT Jacobi |
|---:|---:|---:|---:|
| 1,216 | 28 | 233 | 107 |
| 4,736 | 33 | 592 | 136 |
| 18,688 | 36 | 2000 (failed) | 159 |

At 18,688 unknowns, the new residual is `9.48e-9`, versus `2.11e-8` for
the baseline. The corresponding solve times are approximately 0.61 s and
5.68 s; macro patches take 0.18 s.

With two pre/post steps on the triangle, the new smoother takes 92 iterations
at 2,400 unknowns and 103 at 9,408, compared with 583 and 2000 (failed) for
patches alone. Macro patches take 21 at both sizes. At 9,408 unknowns, the extra
sweeps increase the new smoother's solve time from about 0.27 s to 0.31 s despite
reducing iterations.

### Refinement dependence: the main finding

The updated split smoother still degrades under refinement. At fixed `1e-8`
tolerance, its iteration counts are 37, 83, 110, 134, and 152 from 48 through
9,408 unknowns. The increments become smaller, but these data do not establish
a plateau. The old patch-only counts grow much faster: 41, 141, 362, 856, then
2000 with tolerance missed.

One further triangle refinement gives 37,248 unknowns. The new smoother hits
2000 iterations with true residual `5.68e-8`, failing the requested `1e-8`.
To distinguish this stringent-tolerance failure from the earlier iteration
trend, the following additional runs used the same hierarchy:

```sh
OPENBLAS_NUM_THREADS=1 addisons_experiments/build/mg_compare -r 6 -smoother macro \
  > addisons_experiments/output/mg_macro_deeper.csv
OPENBLAS_NUM_THREADS=1 addisons_experiments/build/mg_compare -r 6 -smoother split -tol 1e-7 \
  > addisons_experiments/output/mg_triangle_relaxed.csv
```

At the fixed, relaxed `1e-7` tolerance, the split-plus-HCT iteration counts are:

| JM unknowns | Split + HCT Jacobi iterations |
|---:|---:|
| 15 | 1 |
| 48 | 37 |
| 168 | 74 |
| 624 | 100 |
| 2,400 | 121 |
| 9,408 | 138 |
| 37,248 | 153 |

The last refinement therefore raises iterations from 138 to 153 (about 11%)
while unknowns increase almost fourfold. Solve time rises from about 0.24 s to
1.32 s (about 5.4 times). This is much milder refinement degradation than
patches alone, but neither constant iterations nor demonstrated linear solve
cost. Do not compare counts across different tolerances as one sequence.

The macro control also fails `1e-8` at 37,248 unknowns, ending with residual
`6.49e21` after 2000 iterations. Thus the strict-tolerance failure at this depth
is not unique to the new HCT smoother. These endpoint measurements do not
identify its cause: numerical accuracy, continuing CG after residual stagnation,
and hierarchy stability require further investigation. In particular, the
failure cannot be attributed solely to inadequate divergence-free smoothing.
The macro command returns 3; the relaxed split command returns 0.

**Conclusion:** adding local divergence-free corrections greatly reduces
refinement degradation and resolves the split-patch failures over the original
mesh range. Iterations still increase with refinement, and an extra level fails
the strict tolerance. Refinement-independent convergence is not established.
Over the original range, the combined smoother still requires substantially
more iterations and time than macro patches. The results support correcting the divergence-free modes lost
by splitting macro patches, but do not prove a stable kernel decomposition.

## Validation

All three standalone CTest entries pass. The focused tests check:

- HCT Jacobi action against an explicit sum over basis-image corrections,
  normalized by their energy under the full JM mass-plus-div-div operator;
- zero divergence of each HCT basis image;
- symmetry and sampled positive quadratic action of all four V-cycle variants;
- macro-patch agreement between bases and exactly unchanged moment transfer.

These finite regression checks supplement the damping bound; they are not a
mesh-uniform convergence proof.
