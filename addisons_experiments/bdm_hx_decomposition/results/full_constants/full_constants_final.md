# Full squared product-norm constants, alpha = 1: final numerical report

Worker model: `claude-opus-5-5`. Compact data:
`results/full_constants/full_constants_summary.json`.

## Definitions (as in FULL_CONSTANTS.md, unchanged)

The setting is the unit square with n x n squares, each cut into two triangles,
and BDM1. The input norm is A(u) = ||u||^2_H(div).

- **CoptHX(h)** = lambda_max(B^{-1}, A). It is the optimal discrete HX constant over
  u_h = s_h + R q_h + sum_v r_v, measured in ||s_h||^2_H1 + ||q_h||^2_H1 +
  sum_v ||r_v||^2_H(div). Here s_h is continuous vector P1 and q_h is continuous
  P2. The r_v lie in the endpoint vertex patch spaces, which include boundary
  vertices, with no essential BCs. Full H1 norms include the mass terms.
- **Cavg(h)** is the explicit averaged construction applied to the continuous
  minimal preimage (s, q) of u_h. It takes s_h = J_h s (area-averaged trianglewise
  P1 projection), q_h = Pi_Q q (canonical commuting transfer), and
  sum_v r_v = Pi_V s - J_h s, split into the vertex patches. The three energies
  are summed and maximized as **one** quadratic form (generalized eigenproblem).
  Cavg is computed without potential mean normalization. Cavg_mean_subtracted,
  which subtracts the mean of q_h, is reported alongside it.
- The continuous minimizer is approximated in tensor cubic splines of relative
  resolution r, meaning an (n·r)-interval spline grid. Every Cavg is therefore
  labelled by r.
- Ccont(h, r) is the spline approximation of the continuous minimal constant.
  The exact value is <= 1 and tends to 1. It is secondary here.

## Final table

| n | dim BDM | CoptHX | Cavg r=8 | Cavg r=16 | Cavg r=32 |
|---|---:|---:|---:|---:|---:|
| 5 | 170 | 0.9971658422 | 1.0703425565 | 1.0694274618 | 1.0693113489 |
| 10 | 640 | 0.9993780123 | 1.1441598924 | 1.1430557141 | 1.1429155096 |
| 15 | 1410 | 0.9997372594 | 1.1599686976 | 1.1588282090 | 1.1586833683 |
| 20 | 2480 | 0.9998557963 | 1.1655105747 | 1.1643559399 | not run |
| **30** | 5520 | **0.9999373014** | **1.1694346894** | **1.1682691788** | not run |

Mean subtraction changes Cavg by at most 4.4e-16 on every row. The continuous potential q has
mean-operator norm 2.5e-6 or less.

At the n=30, r=16 maximizer the energy components are s_h 1.14163,
q_h 2.24e-5 and patches 0.02662. The CoptHX optimal split at its own maximizer
is s 2.17e-4, q 0.99950 and patches 2.23e-4. The worst case for the explicit
construction is dominated by the smooth part J_h s.

## Resolution diagnostics (successive spline refinements)

| n | Cavg(8)-Cavg(16) | Cavg(16)-Cavg(32) | product-map change 8→16 (smooth / potential / patches) | 16→32 |
|---|---:|---:|---|---:|
| 5 | 9.15e-4 | 1.16e-4 | 3.76e-3 (3.69e-3 / 1.6e-4 / 1.34e-3) | 4.84e-4 |
| 10 | 1.104e-3 | 1.40e-4 | 6.84e-3 (6.78e-3 / 1.8e-4 / 2.11e-3) | 8.74e-4 |
| 15 | 1.140e-3 | 1.45e-4 | 1.024e-2 (1.016e-2 / 1.8e-4 / 2.96e-3) | 1.30e-3 |
| 20 | 1.155e-3 | — | 1.375e-2 (1.366e-2 / 1.8e-4 / 3.82e-3) | — |
| 30 | 1.166e-3 | — | **2.091e-2** (2.079e-2 / 1.8e-4 / 5.53e-3) | — |

The product-map change is the operator norm of the change in the full map
u_h -> (J_h s, Pi_Q q, patch terms), measured in the product norm relative to A.
The row "previous maximizer ratio now" evaluates the r=8 maximizer under the
r=16 form. At n=30 it gives 1.1682639, compared with Cavg(16) = 1.1682692.

- At fixed relative r, the Cavg shift from r=8 to r=16 is nearly independent of n
  (about 1.1e-3 to 1.17e-3). Where r=32 exists (n ≤ 15), the ratio of successive
  shifts is 0.127, which is roughly r^-3 behaviour.
- The product-map change at fixed relative r grows about linearly in n
  (roughly 7e-4 · n). The output map is therefore not uniformly resolved at
  fixed r, even though the output constant is stable. At n=30 the 8→16 map
  change is 2.1e-2, and it is almost all in the smooth component.
- The continuous objective is clearly under-resolved:
  - Ccont(30, 8) = 81.1 and Ccont(30, 16) = 41.3, although the exact value is ≤ 1.
  - The continuous-map change 8→16 at n=30 is 6.31.
  - Even at the Cavg maximizer, the continuous pair energy is 1.0758 at r=8,
    1.0376 at r=16 and (n ≤ 15) about 1.019 at r=32. The observed decrease is roughly
    linear in 1/r; these values alone do not identify its limit.
  - Cavg is much less sensitive, but it is computed from this under-resolved
    minimizer.
- None of these successive differences certifies an error against the
  infinite-dimensional map.

## Extrapolation fits (least squares; diagnostics only)

Each fit uses mesh points n in the subset, all taken at the same labelled
resolution. The exceptions are the rows marked "best r", which mix resolutions:
r=32 for n=5, 10, 15 and r=16 for n=20, 30. The 3-point fits of the 3-parameter
models interpolate exactly. The last column is the maximum residual of the
L+a/n^2 fit.

| series | subset (n) | L+a/n^2 | L+a/n^2+b/n^3 | L+a/n+b/n^2 | max resid (1st) |
|---|---|---:|---:|---:|---:|
| CoptHX | all {5,10,15,20,30} | 1.000048 | 1.000000 | 0.999938 | 4.8e-5 |
| CoptHX | {10,15,20,30} | 1.000012 | 0.999999 | 0.999979 | 5.8e-6 |
| CoptHX | {15,20,30} | 1.000005 | 0.999999 | 0.999989 | 1.2e-6 |
| Cavg r=16 | all | 1.170137 | 1.171838 | 1.173932 | 1.8e-3 |
| Cavg r=16 | {10,15,20,30} | 1.171436 | 1.171412 | 1.171367 | 1.4e-5 |
| Cavg r=16 | {15,20,30} | 1.171425 | 1.171366 | 1.171252 | 1.4e-5 |
| Cavg r=8 | all | 1.171307 | 1.173011 | 1.175108 | 1.8e-3 |
| Cavg r=8 | {10,15,20,30} | 1.172608 | 1.172587 | 1.172546 | 1.4e-5 |
| Cavg r=8 | {15,20,30} | 1.172598 | 1.172541 | 1.172430 | 1.3e-5 |
| Cavg best r | all | 1.170079 | 1.171867 | 1.174085 | 1.9e-3 |
| Cavg best r | {10,15,20,30} | 1.171440 | 1.171517 | 1.171595 | 7.1e-5 |
| Cavg best r | {15,20,30} | 1.171499 | 1.171265 | 1.170818 | 5.3e-5 |

Summary of the fits:

- **CoptHX:** every fit gives 0.99994 to 1.00005, and the larger-mesh fits agree
  with 1 to within 2e-5. This is consistent with the proven lim inf ≥ 1.
- **Cavg, paired r=16 series:** the larger-mesh fits give 1.17125 to 1.17144.
  The all-point fits spread from 1.1701 to 1.1739, and n=5 visibly misfits:
  L+a/n^2 has residual 1.8e-3 there.
- **Cavg, r=8 series:** its intercepts are uniformly about 1.17e-3 higher, which
  matches the per-mesh r-shift.
- **Heuristic r→∞ shift:** if the observed r^-3 ratio (measured only at n ≤ 15)
  persisted, the intercept would move down by a further roughly 1.7e-4. That
  would give about 1.1712 to 1.1713, which is unverified at n ≥ 20.
- **Best-r series:** its {15,20,30} fits are the least trustworthy, because
  mixing r=32 and r=16 injects a spurious 1.4e-4 step between n=15 and n=20.
- These spreads describe variation between models and subsets. They are not
  confidence intervals or enclosures.

## Checks (all rows)

- **CoptHX:**
  - Rayleigh quotient and optimal-split energy agree with lambda_max to 1e-13.
  - Reconstruction of u_h from the optimal split is exact to 1.6e-15 relative.
  - The generalized-eigen residual is 3.4e-12 or less.
  - B^{-1} is exactly symmetric.
- **Cavg (per mesh and r):**
  - The spline stationarity residual is 1.4e-12 or less.
  - The energy identity for the block form has error 1.9e-11 or less at n=30
    (1.5e-10 at n=15, r=32).
  - The reconstruction J_h s + R Pi_Q q + sum r_v = u_h holds to a maximum error
    of 6.0e-9 at n=30, growing mildly with n from 2e-12 at n=5.
  - The Cavg eigen residual is 2e-15 or less, and the Rayleigh quotient agrees to
    1e-15.
- **Dominance:** the minimum eigenvalue of F_avg - B^{-1} is positive in every
  case (1.2e-9 at n=30, r=16), so CoptHX ≤ Cavg holds numerically as a form
  inequality. The scalar comparison CoptHX ≤ Cavg is true everywhere.
- **Potential mean:** mean subtraction changes Cavg by 4.4e-16 or less. The
  continuous q has mean-operator norm 2.5e-6 or less.
- **Continuous form:** the Ccont form decreases monotonically under refinement
  (the minimum change is -3e-9, i.e. rounding-level), as nested minimization
  requires.

## Limitations

- The continuous minimal-preimage map is under-resolved at every r used:
  Ccont is between 2.8 and 81, far above the exact bound of 1. Cavg
  stabilizes much faster than the map, but no computed quantity bounds the error
  of Cavg against the infinite-dimensional minimizer.
- r=32 exists only for n ≤ 15. Only r=8 and r=16 are available at every mesh, so
  the paired r=16 series is the primary Cavg sequence.
- **Feasibility limit at n=30:**
  - The r=16 run peaked at 11.4 GB RSS (process maximum; 10.9 GB recorded
    in-driver). The spline dimension was 233,289, with 3.9e8 LU entries, and the
    run took 47.8 min of wall time.
  - n=30 r=32 would need 923,521 spline dofs and an estimated 1.7e9 or more LU
    entries (above 20 GB). That is infeasible on this 22 GB machine, which had
    about 13 GB free.
  - n=20 r=32 and n=40 r=16 (413k dofs, dense N up to 9760) are estimated at
    roughly 10 GB of LU plus dense eigen storage. They were not attempted, per
    the no-unbounded-sweeps instruction.
  - The nearest feasible paired meshes are therefore n=30 at r=8 and r=16, which
    were completed as requested.
- The fits use 3 to 5 points with n ≤ 30. Their intercepts are extrapolation
  diagnostics only.
- The n ≤ 20 results were reused from the earlier runs of the same driver; they
  were not recomputed.

## Implementation changes made in this wrap-up

1. `full_constants.py` received memory-only edits. The original is saved at
   `agent_artifacts/full_constants_prewrapup.py.bak`.
   - (a) The LU fill is recorded as `lu_storage_nnz = lu.nnz` instead of
     `lu.L.nnz + lu.U.nnz`, because the latter copies the factor.
   - (b) All spline-LU-dependent quantities (the Ccont form, mean row and
     direct evaluations) are computed first. The spline LU is then freed before
     the remaining dense eigensolves. The output record keys and their order are
     otherwise unchanged.
   - (c) Peak RSS is recorded (`peak_rss_gb`).
   - Verification: rerunning n=5 at r = 8, 16, 32 gave every numeric field
     bitwise identical to the prior JSON, apart from the renamed LU count.
2. New script `fit_full_constants.py`. It reads `agent_artifacts/full_constants_n*.json`
   and writes `results/full_constants/full_constants_summary.json`, which holds
   the compact per-mesh data, series, fits (with the exact points used) and
   provenance.

No changes were made to the MFEM library, to the parent-authored markdown, or
to git.

## Reproduction

```
cd addisons_experiments/bdm_hx_decomposition
python3 full_constants.py --mesh 5  --refinements 8 16 32
python3 full_constants.py --mesh 10 --refinements 8 16 32
python3 full_constants.py --mesh 15 --refinements 8 16 32
python3 full_constants.py --mesh 20 --refinements 8 16
systemd-run --user --scope -p MemoryMax=15G -p MemorySwapMax=1G \
  python3 full_constants.py --mesh 30 --refinements 8 16 --chunk 128
python3 fit_full_constants.py
```

Outputs go to `agent_artifacts/full_constants_n{N}.json`, together with maximizer
.npy files and logs, all ignored. `--chunk` only changes the column blocking of the spline solves (default 256), and hence the diagonal blocks used by the energy-identity check.

Provenance:

- Model: claude-opus-5-5.
- git HEAD: 8c56dfdaef (branch johnson-mercier). The experiment scripts are
  uncommitted.
- Software: Python 3.14.7, numpy 2.4.6, scipy 1.16.2 (SuperLU), FlexiBLAS.
- Machine: AMD Ryzen 7 5825U, 22 GB RAM, Linux 7.1.13-200.fc44.
- Run date: 2026-09-25.
