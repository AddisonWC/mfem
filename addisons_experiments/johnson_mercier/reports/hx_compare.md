# Johnson–Mercier HX comparisons

These are historical measurements from before the directory reorganization.
Program names and commands below use the current standalone layout; see the
[README](../../README.md) for build instructions and provenance.


`hx_compare` solves the symmetric-tensor H(div) problem

\[
 a(\sigma,\tau)=(\sigma,\tau)+(\operatorname{div}\sigma,
                                      \operatorname{div}\tau).
\]

It compares two JM coordinate bases, three smoothers, and two auxiliary H1
spaces. The HCT Airy correction is present in every configuration. There are no
essential stress boundary conditions. The current example requires a connected,
conforming, straight 2D triangle mesh. It is a serial, assembled experiment.

## Build and run

Build the standalone project following [the README](../../README.md). SuiteSparse is recommended:
auxiliary inverses use UMFPACK when available. Without it, a fixed dense LU
factorization is used, limited to auxiliary matrices of size 2000. There are no
tolerance-based inner iterative solves.

From `addisons_experiments/build`:

```sh
OPENBLAS_NUM_THREADS=1 ./hx_compare -r 2 -levels 4 -repeat 5 > triangle.csv
OPENBLAS_NUM_THREADS=1 ./hx_compare -m data/inline-tri.mesh -r 0 -levels 3 > square.csv
./hx_compare -r 3 -basis vertex -smoother jacobi -h1 split
```

The selectors accept:

| Option | Choices |
|---|---|
| `-basis` | `all`, `moment`, `vertex` |
| `-smoother` | `all`, `jacobi`, `macro`, `split` |
| `-h1` | `all`, `unsplit`, `split` |

`-r` is the initial refinement count; `-levels` runs that level and successive
uniform refinements. `-repeat` controls repeated solves for median timings.
`-damping` weights only the smoother correction (default 1), relative to the two
auxiliary corrections. For comparison with the older example's patch scaling,
use `-damping 0.33`. A seeded random load in the common moment coordinates is
the default; `-smooth-rhs` instead assembles a nonconstant smooth tensor load.

Standard output is CSV, and options/backend information go to standard error.
Every solve starts at zero and uses the same load and original moment-coordinate
operator. CG stops using a common unpreconditioned relative residual, with the
true residual recomputed before accepting convergence. The exit status is 3 if
any configuration fails to converge. `residual` reports the final recomputed
relative residual.

## Basis and patch meaning

The existing `JohnsonMercierFECollection()` remains unchanged by default. The
new selection is:

```cpp
JohnsonMercierFECollection fec(JMBasis::SplitVertex);
```

Its serialized collection name is `JM_2D_P1_SplitVertex`. It has the same space
and entity ownership as `JM_2D_P1`: four DOFs per edge, three per triangle, zero
vertex-owned DOFs. The edge coordinates are `(nn,nt)` traction values at each
endpoint in the edge's normal/tangent frame. The interior coordinates are the
physical tensor components `(s00,s01,s11)` at the barycenter. Reversing an edge
swaps its two endpoint pairs without signs.

The dual basis splits into local vertex subspaces of dimensions 4,4,4,3. The
interior basis functions are the barycenter hat times a constant symmetric
tensor. The finite element's interpolation remains the original moment
interpolation, expressed in the selected coordinates; it is not replaced by
point interpolation of arbitrary source functions. Refinement transfer retains
the same moment interpolant, including the existing physical mapping correction.
JM spaces on uniformly refined macro meshes are not generally nested.

`GetMomentToSplitVertexMatrix` exposes the local physical coefficient change.
The example assembles its inverse globally as `B`, so moment coefficients equal
`B * vertex_coefficients`. For a vertex-basis inverse `R`, its action in the
common moment coordinates is `B R B^t`. This keeps the operator, right-hand side,
and convergence criterion identical in all comparisons.

The smoothers are:

- `jacobi`: the inverse diagonal in the selected basis.
- `macro`: an exact solve on each original vertex patch, collecting all four
  DOFs of incident edges and all three DOFs of incident triangles. These patches
  overlap. Their subspaces are the same in both bases; the implementation uses
  the selected coordinates and verifies their equivalence in a unit test.
- `split`: exact solves on the intrinsic split-vertex spaces. An original vertex
  gets two DOFs per incident macro edge; each barycenter gets a separate block of
  three. In vertex coordinates, these blocks partition the DOFs. In moment
  coordinates the smoother is applied through `B`, not by selecting moment DOF
  indices. Consequently its two basis rows intentionally describe the same
  preconditioner.

Patch matrices are factored once, without explicitly forming their inverses.
Scalar Jacobi depends on the basis directions; exact patch solvers depend only
on their subspaces. The weighted-mass condition-number bound of 4 applies to
exact split-vertex patch solves for that mass operator. It does not establish a
bound of 4 for this div-div-plus-mass HX experiment. General trace/deviatoric
compliance weights are not implemented by this example.

## Split auxiliary H1 space

The optional explicit Alfeld mesh retains the original vertex numbers, appends
one barycenter per macro triangle, and creates three child triangles. Original
boundary segments and attributes are preserved. Ordinary vector H1 assembly
uses Frobenius component weights `(1,2,1)` for both mass and diffusion.

Continuous split P1 tensors are already JM functions. Their inclusion is built
from edge endpoint traction values and barycenter tensor values, then mapped to
moment coordinates. This dedicated inclusion avoids pairing unrelated element
indices through `DiscreteLinearOperator`. Unit tests check both the actual
piecewise fields and the identity between the split H1 mass and the restricted
JM mass.

The split H1 space adds three DOFs per macro triangle to the unsplit space. Its
additional directions are precisely the barycenter patch spaces, but they are
coupled through the auxiliary H1 solve. A larger auxiliary space can therefore
help Jacobi without necessarily benefiting the exact patch smoothers as much.

The HCT space stays on the macro mesh. Its affine kernel is fixed by setting the
value and both first derivatives at one vertex; the corresponding residual
entries are also zeroed before applying the inverse.

## Timing interpretation

The CSV separates moment assembly, vertex assembly, basis-map construction,
auxiliary setup, smoother setup, median solve time, and median time for one HX
application. All times are wall-clock seconds. Auxiliary setup is shared across
the six smoother/basis rows for each H1 choice; its time is repeated in those
rows and must not be summed as if setup were performed six times.

The harness assembles both JM bases even for a single selected configuration.
Those costs are reported separately; both assemblies would not be necessary in
a production solver using only one basis. Vertex-basis timings currently include
the explicit sparse coordinate maps during application. Native matrix assembly
also evaluates the physical basis correction at quadrature points, so its cost
is relevant separately from iteration savings.

Use an idle machine, a fixed BLAS thread count, and several repetitions. Small
timing differences, particularly between identical split-patch rows, measure
timing variability rather than different algorithms. Direct auxiliary solves
isolate these decomposition comparisons; the example does not establish
scalability of the auxiliary inverses themselves.

## Recorded experiments

[hx_results.csv](../results/hx_results.csv) records 144 configurations, each with
five repeated solves. All converged with a recomputed relative residual below
`1e-8`. These measurements used a Release build with GCC 16.2.1, UMFPACK, and
`OPENBLAS_NUM_THREADS=1` on an AMD Ryzen 7 5825U. The dataset names correspond to:

```sh
# Run in addisons_experiments/build with OPENBLAS_NUM_THREADS=1.
./hx_compare -r 2 -levels 4 -repeat 5
./hx_compare -r 2 -levels 4 -repeat 5 -smooth-rhs
./hx_compare -m data/inline-tri.mesh -r 0 -levels 3 -repeat 5
./hx_compare -r 5 -repeat 5 -damping 0.33
```

They are labeled `triangle_random`, `triangle_smooth`, `square_random`, and
`triangle_damped`. The random triangle dataset was repeated to investigate
timing variation; the CSV retains the final run.

With the random RHS and damping 1, the finest-grid iteration counts were:

| Smoother / basis | Triangle: unsplit H1 | Triangle: split H1 | Square: unsplit H1 | Square: split H1 |
|---|---:|---:|---:|---:|
| Jacobi / moment | 117 | 73 | 103 | 68 |
| Jacobi / vertex | 52 | 49 | 47 | 45 |
| Macro patches / either | 28 | 27 | 26 | 25 |
| Split patches / either | 36 | 36 | 34 | 34 |

The triangle has 1,024 macro elements and 9,408 JM DOFs; its H1 spaces have
1,683 and 4,755 DOFs. The square has 512 macro elements and 4,736 JM DOFs; its
H1 spaces have 867 and 2,403 DOFs. Exact patch counts agreed between the two
bases at every tested level.

For the triangle random RHS with unsplit H1, refinement gave:

| JM DOFs | Moment Jacobi | Vertex Jacobi | Macro patches | Split patches |
|---:|---:|---:|---:|---:|
| 168 | 40 | 34 | 19 | 24 |
| 624 | 66 | 42 | 21 | 29 |
| 2,400 | 93 | 49 | 24 | 34 |
| 9,408 | 117 | 52 | 28 | 36 |

The new basis substantially improves Jacobi on the finer tested meshes. Exact
macro patches still need the fewest iterations, and split patches lie between
macro patches and vertex Jacobi. These finite-level measurements do not prove
mesh-independent iteration bounds.

The smooth triangle load gives the same overall picture at 9,408 DOFs:
moment/vertex Jacobi take 99/41 iterations with unsplit H1 and 62/40 with split
H1. Macro patches take 23 versus 26, and split patches take 30 in both cases.
The larger H1 space primarily helps moment Jacobi in these experiments; it is
not uniformly better once vertex coordinates or exact patches are used.

For the finest square problem, median solve times in milliseconds were:

| Configuration | Unsplit H1 | Split H1 |
|---|---:|---:|
| Moment Jacobi | 80.3 | 61.8 |
| Vertex Jacobi | 39.3 | 43.2 |
| Macro patches, moment coordinates | 27.4 | 30.0 |
| Split patches, vertex coordinates | 30.1 | 35.3 |

Split-patch setup was about 0.57 ms versus 3.34 ms for macro patches with
unsplit H1. Their cheaper setup and application partly compensate for the
additional iterations. Native vertex-basis matrix assembly took about 71 ms
versus 13 ms for moment assembly on this mesh: optimizing or precomputing the
physical basis correction remains useful for applications where setup dominates.

The triangle random-run timings exhibited substantial variation even between
the mathematically identical split-patch rows (roughly 95 versus 163 ms on the
finest unsplit-H1 case). This persisted on repeating that dataset, so those
timings should not be used to rank the equivalent representations. The raw
measurements are retained; the iteration counts are the stronger comparison.

Damping is a separate parameter. At 9,408 DOFs, reducing it from 1 to 0.33
changes macro-patch counts from 28 to 24 with unsplit H1 and 27 to 23 with split
H1. It changes split-patch counts from 36 to 44 and 36 to 38, respectively.
The old example's factor 0.33 therefore should not be copied indiscriminately
to the new patches.

For subsequent experiments, vertex Jacobi with unsplit H1 is a useful inexpensive
baseline. Compare it with both patch choices, measure setup separately, and
tune damping per smoother before selecting a default. On these problems, exact
macro patches remain the strongest option by iteration count and solve time;
the split patches offer substantially smaller local factorizations.

The original validation recorded 3,215 passing assertions in the JM/HCT/HX tests, covering
physical basis equivalence, vertex support, traction continuity, canonical
interpolation, refinement/update transfer, Airy divergence, split H1 inclusion,
macro-patch invariance, and mass decomposition bounds. The existing `ex43` and
`ex43_hx` CTest cases and the comparison smoke test also passed. The dense inverse
fallback was separately exercised on all 12 small-mesh configurations.
