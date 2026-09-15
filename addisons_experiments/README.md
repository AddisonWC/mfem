# Addison's experiments

Research drivers, notes, and selected results for this development branch.
Useful pieces may eventually become library features or maintained miniapps;
the remainder can be discarded by deleting this directory. MFEM's build and
normal tests do not reference it.

## Build and test

First build this branch of MFEM with CMake. From the repository root, using the
existing `build-hx` configuration:

```sh
cmake --build build-hx --target mfem -j 4
cmake -S addisons_experiments -B addisons_experiments/build \
  -DMFEM_DIR="$PWD/build-hx"
cmake --build addisons_experiments/build -j 4
ctest --test-dir addisons_experiments/build --output-on-failure
```

`MFEM_DIR` can also name an installed MFEM CMake package built from this branch.
The experiments require its split-vertex JM support. `MFEM_SOURCE_DIR` defaults
to the parent directory and supplies mesh inputs and Catch/unit-test headers.
No external test dependency is downloaded. Set `BUILD_TESTING=OFF` to build
only the drivers. Rebuild MFEM separately after changing library code.

SuiteSparse is recommended for larger runs. Without it, the experiments use
their size-limited dense direct solver. The default mesh is copied into the
experiment build and its absolute path compiled into each driver; explicit
`-m` paths are relative to the invocation directory.

```sh
cd addisons_experiments/build
OPENBLAS_NUM_THREADS=1 ./hx_compare -r 2 -levels 4 -repeat 5 > hx.csv
OPENBLAS_NUM_THREADS=1 ./mg_compare -r 3 > mg.csv
```

## Contents

- `johnson_mercier/common.hpp`: shared experimental meshes, maps, and solvers.
- `johnson_mercier/{hx,mg}_compare.cpp`: comparison drivers.
- `johnson_mercier/tests/`: experiment validation, run only by this project.
- [HX report](johnson_mercier/reports/hx_compare.md) and
  [multigrid report](johnson_mercier/reports/mg_compare.md).
- `johnson_mercier/results/`: selected CSV snapshots, including the updated HCT-Jacobi multigrid runs.

The original reports and CSVs were moved from `examples/`; the historical CSVs
remain unchanged. The multigrid report now also records the split-patch plus
HCT-Jacobi smoother, with new measurements in `mg_hct_results.csv`.
HX materials came from commit `3f34c7b398`, and multigrid materials from
`81fce13057`. These identify the commits that recorded the results, not necessarily
the exact revisions used to run them. The HX report retains its historical measurements and limitations. See the
multigrid report for reproduction commands and provenance of the new runs.

For future snapshots, record the source revision (and uncommitted changes), exact
commands, compiler/build settings, solver backend, hardware, and thread count.
Keep routine outputs in ignored `build/` or `output/` directories. Preserve
library-level basis, interpolation, orientation, and transfer tests in MFEM's
normal test suite, independently of this project's solver tests.
