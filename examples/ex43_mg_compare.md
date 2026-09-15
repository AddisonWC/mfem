# Johnson–Mercier multigrid with split patches

`ex43_mg_compare` investigates geometric multigrid for the same stress problem as
`ex43` and `ex43_hx_compare`:

\[
 A(\sigma,\tau)=(\sigma,\tau)+(\mathrm{div}\,\sigma,\mathrm{div}\,\tau).
\]

There are no essential stress boundary conditions. The example is serial and
assembled, on a connected, straight, conforming triangular mesh. It uses a
symmetric V-cycle with equal pre/post smoothing and a fixed direct coarse solve.
It reports each prefix of the hierarchy, including the coarse-only solve.

## Prolongation is fixed

The hierarchy is constructed **once using the original moment JM collection**.
Every smoother comparison borrows the identical prolongation objects from this
hierarchy. Restriction is their transpose. Level operators are rediscretized in
moment coordinates, as in `ex43`; they are shared across all comparisons.
JM macro refinements are generally nonnested, so this is the existing moment
interpolation transfer, not a claim of exact inclusion or Galerkin coarsening.

Only the smoother uses the split basis. If `B` maps split coefficients to moments,
its matrix is `Av = B^T A B`, and the split smoother acts as

\[
 S_m r = \omega B\left[\sum_p I_p
 (I_p^T A_v I_p)^{-1}I_p^T\right]B^T r.
\]

The coarse correction still uses the original `P`, without inserting `B` into
prolongation. This is essential: selecting the split patch DOF indices directly
in the moment matrix would describe different subspaces.

The three configurations are:

| Name | Patch subspaces |
|---|---|
| `macro` | Original vertex stars, in moment coordinates |
| `macro-vertex` | The same vertex stars, factored in split coordinates; a control for basis invariance |
| `split` | Two traction DOFs per incident edge endpoint, plus separate three-DOF barycenter blocks |

The macro configurations should agree to roundoff. Split patches are smaller,
but a basis change alone does not guarantee that their smoothing controls the
error left by the coarse correction. In particular, the mass-only patch bound
from the HX investigation does not give a multigrid bound for div-div plus mass.

Patch corrections are additive and damped by `-damping` (default 0.33).
The example restricts this to `(0,0.5)`: at most four patch subspaces contribute
on an element (three for macro patches), so the elementwise Cauchy–Schwarz bound
makes this a conservative common interval for positive symmetric smoothing.

## Build and run

```sh
cmake --build build-hx --target ex43_mg_compare -j 4
cd build-hx/examples
OPENBLAS_NUM_THREADS=1 ./ex43_mg_compare -r 5 > triangle.csv
OPENBLAS_NUM_THREADS=1 ./ex43_mg_compare -m ../data/inline-tri.mesh -r 3 > square.csv
OPENBLAS_NUM_THREADS=1 ./ex43_mg_compare -r 5 -s 2 > two-sweeps.csv
```

The configured examples Makefile also supports `make ex43_mg_compare`.
`-cr` sets refinements before the coarse grid is constructed; `-r` adds levels
above it. `-smoother` selects one configuration or `all`. `-s` controls the number
of smoothing steps on each side of the coarse correction. SuiteSparse supplies
UMFPACK when enabled; otherwise the shared exact solver uses dense LU and limits
the coarse matrix to 2000 unknowns.

All solves start from zero with the same seeded random moment-coordinate load
at a given depth. PCG uses a common unpreconditioned relative residual tolerance
(default `1e-8`), verified by recomputing `Ax-b`. Failure returns exit status 3.
CSV `setup_s` measures smoother factorization and V-cycle construction only;
shared mesh, transfer, operator, coordinate-map and coarse-factorization costs
are excluded. `solve_s` is one solve, so small timing differences are not reliable.

## Interpretation

Unlike HX, this V-cycle has no separate HCT Airy correction. Div-div has a large
divergence-free kernel. A useful smoother must correct fine-scale modes in or
near that kernel that the coarse space misses. Splitting a macro patch into
small endpoint and barycenter blocks removes coupled local corrections.
Rapid iteration growth with these blocks is consistent with poor smoothing of
those modes; the experiment alone does not prove a kernel-decomposition theorem.

The validation case in `tests/unit/fem/test_jm_hx.cpp` checks that the full
macro-patch V-cycle agrees between bases, each of the three V-cycles is symmetric
and has positive quadratic action on a test vector, the split cycle differs
from the macro cycle, and the shared moment prolongation action remains exactly
unchanged. These sampled symmetry/positivity checks are regression checks, not
spectral bounds.

## Recorded results

[ex43_mg_results.csv](ex43_mg_results.csv) contains 63 configurations from a
Release build with UMFPACK and `OPENBLAS_NUM_THREADS=1`. The first three datasets
use the commands above; `triangle_damped` uses `-r 4 -damping 0.49`. The tolerance
is `1e-8` and the iteration limit is 2000 throughout. Failed solves are retained
with `converged=0`; the first three commands therefore return exit status 3.

On the reference triangle with one pre/post step and weight 0.33:

| JM unknowns | Macro patches (either basis) | Split patches |
|---:|---:|---:|
| 15 | 1 | 1 |
| 48 | 13 | 41 |
| 168 | 20 | 141 |
| 624 | 27 | 362 |
| 2,400 | 32 | 856 |
| 9,408 | 34 | 2000 (failed tolerance) |

At 9,408 unknowns the split run's recomputed relative residual is `2.12e-8`.
This is a failure of the requested stopping criterion, not a claim that the
iterate diverges. The large growth in iterations is already visible at the
smaller levels that do meet tolerance. Macro and macro-vertex iteration counts
agree at every recorded level and setting.

On the square (`inline-tri.mesh`), the nontrivial levels have 1,216, 4,736 and
18,688 unknowns. Macro patches take 28, 33 and 36 iterations; split patches take
233, 592 and 2000, with the last run missing tolerance (`2.11e-8` residual).

Increasing the smoothing count to two steps on each side reduces the triangle
counts at 2,400 unknowns from 32 to 21 for macro patches and from 856 to 583 for
split patches. At 9,408 unknowns macro patches take 21 iterations, while split
patches again reach the limit (`2.07e-8` true residual). Increasing the weight to
0.49 gives 33 versus 689 iterations at 2,400 unknowns. These changes improve some
split results but do not remove their strong refinement dependence.

**Conclusion:** for this unchanged-transfer geometric V-cycle, the split patches
are substantially weaker smoothers than the original macro patches. Their
smaller factorizations do not compensate for the iteration growth in these
experiments. The favorable split-patch behavior in the HX experiment does not
carry over to this multigrid coarse correction.

The new example's CTest smoke test passes, as do 415 assertions in the two
focused JM HX/multigrid unit cases, including the unchanged-transfer regression.
