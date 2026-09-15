# Working on experiments

- Put research drivers, solver comparisons, reports, and selected results here.
- Keep this project standalone: MFEM source, build files, and normal tests must
  never depend on this directory. Deleting it must require no other edits.
- Keep reusable finite-element capabilities in MFEM's existing library layout,
  with independent regression tests in tests/unit/fem.
- Build and test experiments explicitly through this directory's CMake project.
- Put routine generated output in build/ or output/ (both ignored). Commit only
  selected results with reproduction commands and available provenance.
- When practical, commit reusable library changes separately from experiments.
