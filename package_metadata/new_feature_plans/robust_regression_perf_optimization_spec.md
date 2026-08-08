# Robust-Regression Bootstrap-Loop Perf Optimization Spec

Generated: 2026-07-29

## Scope

This spec defines a profile-first performance optimization pass on
`compute_robust_rand_bootstrap_parallel_cpp` in
[../EDI/src/fast_robust_regression.cpp](../EDI/src/fast_robust_regression.cpp),
the OpenMP-parallel loop that calls `fast_robust_regression_internal` once per
randomization/bootstrap replicate for `InferenceContinRobustRegr`'s smoothed-BRT
path (`compute_fast_rand_bootstrap_distr`,
[../EDI/R/inference_continuous_robust_regr.R:330](../EDI/R/inference_continuous_robust_regr.R)).

As a secondary, lower-priority check, this spec also covers the single-shot
`fast_robust_regression_cpp` / `fast_robust_regression_weighted_cpp` entry
points, now used directly (not through the parallel bootstrap kernel) by
`InferenceContinKKRobustRegrOneLik` and `InferenceContinKKRobustRegrIVWC`
after the 2026-07-29 `use_rcpp` wiring — this pass must not regress those.

Related documents:

- [perf_experiments_final.md](perf_experiments_final.md) — do not duplicate
  perf work already tracked there.
- [../package_tests/path_audits_source.R](../package_tests/path_audits_source.R)
  line 34 — `InferenceContinRobustRegr` audit row: `skip_brt="ci"`, i.e. the
  BRT confidence-interval path is still considered too slow to run in the
  comprehensive test suite even after the existing smoothed-noise fast-path
  optimization.

## Motivation

`compute_robust_rand_bootstrap_parallel_cpp` calls
`fast_robust_regression_internal` fresh, from scratch, once per replicate
(B is typically 251–501 for BRT paths in this codebase), with
`smart_cold_start = true` and no state carried across replicates beyond what
each thread's stack provides. For the small `p` (2–5 covariates) typical of
KK matched-pair / simple robust-regression designs, per-call heap allocation
(`b_free`, `r`, `abs_r`, `XtWX`, a fresh `Eigen::LDLT` workspace, etc. —
all allocated new on every single replicate call) plausibly dominates wall
time far more than any remaining un-vectorized arithmetic.

This matters because the GEMM/crossprod core has **already** had a real SIMD
pass: `weighted_crossprod` / `weighted_crossprod_rhs` in
[../EDI/src/_helper_functions.h](../EDI/src/_helper_functions.h) already use
explicit `#pragma omp simd` inner loops, and `symmetric_crossprod` already
uses BLAS `DSYRK`. Blindly hand-vectorizing further in
`fast_robust_regression.cpp` without profiling risks spending effort on code
that is not the bottleneck. `path_audits_source.R:34`'s `skip_brt="ci"` is
concrete evidence that a real perf problem remains somewhere in this path —
this spec's job is to find out where before deciding how to fix it.

## Non-Goals

- Do not modify `weighted_crossprod`, `weighted_crossprod_rhs`, or
  `symmetric_crossprod` in `_helper_functions.h` unless the profiling phase
  (Section "Profiling Phase") specifically implicates them as the dominant
  cost. They already have a dedicated SIMD/BLAS optimization pass.
- Do not change any statistical behavior, estimate, or standard error this
  kernel produces. This is a wall-clock-only pass; every optimization must be
  validated as numerically identical (within floating-point tolerance) to
  the pre-optimization output on a fixed seed.
- Do not commit to SIMD intrinsics as the deliverable. If profiling shows
  allocation/dispatch overhead dominates, the deliverable is scratch-buffer
  reuse, not hand-written AVX code — see "Candidate Optimization
  Directions."
- Do not touch the `use_rcpp` wiring in `InferenceContinKKRobustRegrOneLik`,
  `InferenceContinKKRobustRegrIVWC`, or `InferenceContinRobustRegr`
  (2026-07-29 change) beyond verifying it doesn't regress.

## Profiling Phase

Use this repo's existing `EDI_DEBUG_SYMBOLS=1` build flag (`Makevars`,
already wired for `perf annotate`-style profiling: keeps `-O3` but adds
`-g1 -fno-omit-frame-pointer`) to build a debug-symboled `.so`, then profile
a representative call to `compute_fast_rand_bootstrap_distr` (B = 501,
realistic n and p for a KK matched-pair robust-regression design) with
`perf record` / `perf report` (or equivalent) to get a real breakdown of
where time goes inside a single replicate:

1. Cold-start QR solve (`Eigen::ColPivHouseholderQR<...>::solve`)
2. MAD scale estimation (`std::sort` / `std::nth_element` on `abs_r`)
3. Per-IRLS-iteration weight computation (Huber/bisquare `.select()` array
   ops)
4. `weighted_crossprod(X_free, res.w)` / `weighted_crossprod_rhs(...)`
5. `Eigen::LDLT<Eigen::MatrixXd>` construction + solve
6. Everything else (heap allocation, `RobustModelResult` construction,
   OpenMP dispatch overhead)

Report the percentage of wall time attributable to each bucket. This
breakdown determines which of the candidate directions below (if any) are
worth implementing — no direction should be implemented before the profile
confirms it is a real contributor.

## Candidate Optimization Directions

These are candidates, not commitments. Implement only what the profile
justifies, roughly in this priority order based on the reasoning in
"Motivation":

1. **Per-thread scratch-buffer reuse across the `nsim` loop** (most likely
   winner). Currently `fast_robust_regression_internal` allocates `b_free`,
   `r`, `abs_r`, per-iteration `XtWX`/`XtWy`, and a fresh `LDLT` workspace on
   every call. `compute_robust_rand_bootstrap_parallel_cpp` already hoists
   `y_b`/`X_b` outside the `for (b ...)` loop per OpenMP thread — the same
   pattern should extend to `fast_robust_regression_internal`'s internals,
   e.g. via an optional pre-allocated workspace struct passed by pointer,
   reused across all replicates a thread processes.
2. **Small-`p` specialized solve.** If `LDLT` construction/solve is a real
   contributor at `p` = 2–5, consider a templated or hand-written small-matrix
   solve path instead of the fully generic `Eigen::LDLT`.
3. **Avoid full sort for MAD when possible.** `std::sort`/`std::nth_element`
   run once per replicate (not per IRLS iteration). Only pursue this if the
   profile shows it as a meaningful fraction of per-replicate time; a
   partial-selection approach (`nth_element` is already used for n ≥ 512,
   full sort below that) may already be near-optimal for the small `n`
   typical here.
4. **Explicit SIMD intrinsics.** Lowest priority — only for a scalar loop the
   profile identifies as hot and that the compiler is demonstrably failing
   to auto-vectorize (check `-march=native` codegen / `-fopt-info-vec` before
   hand-writing intrinsics).

## `estimate_only` Path Handling

The existing code already cleanly separates the two paths:
`fast_robust_regression_cpp`'s `estimate_only = TRUE` branch returns
immediately after the IRLS loop, before `ssq_j`/`fisher_information` are
computed (`fast_robust_regression.cpp`, the early `if (estimate_only) {
return ... }` block). `compute_robust_rand_bootstrap_parallel_cpp` always
calls with `estimate_only = true`. Any scratch-buffer or allocation-reuse
work must preserve this split exactly: buffers sized for the SE/vcov path
(`XtWX` expansion, `ssq_j` computation) must not be allocated or computed
when `estimate_only = TRUE`, matching current behavior.

## Testing

- Regression baseline: `EDI/tests/testthat/test-brt-smoothed-robust-regr-kernel.R`,
  `test-robust-regression-mad.R`, `test-robust-regression-xtx-cache.R`, and
  the 2026-07-29 `test-kk-robust-regr-use-rcpp.R` must all continue passing
  unchanged.
- Add a fixed-seed before/after numerical-parity test:
  `compute_robust_rand_bootstrap_parallel_cpp` output on a fixed
  `(y0, Xc, i_mat, w_mat, delta, method, noise_mat)` input must match the
  pre-optimization output within `1e-10` (this is a refactor of allocation
  strategy, not algorithm, so parity should be near machine-precision, not a
  loose statistical tolerance).

## Benchmark

Profile before (see "Profiling Phase"), implement the justified direction(s),
profile after. Report:

- Wall-clock for `compute_fast_rand_bootstrap_distr` at B = 251 and B = 501,
  at a KK matched-pair n/p representative of `path_audits_source.R`'s
  `InferenceContinRobustRegr` audit row.
- Wall-clock for the standalone `fast_robust_regression_cpp` /
  `_weighted_cpp` single-call path (non-regression check for the KK classes'
  `use_rcpp` wiring).
- Whether `skip_brt="ci"` in `path_audits_source.R:34` can be lifted as a
  result of this work — this is the concrete success signal the motivation
  section is chasing, not a specific speedup multiplier.

## Implementation Phases

### Phase 1: Profiling
Build with `EDI_DEBUG_SYMBOLS=1`, profile `compute_fast_rand_bootstrap_distr`
at realistic B/n/p, produce the per-bucket breakdown described above.

### Phase 2: Targeted optimization
Implement only the candidate direction(s) the profile justifies. Expect this
to be scratch-buffer reuse (Candidate 1) based on the reasoning in
Motivation, but the profile is authoritative, not this expectation.

### Phase 3: Correctness verification
Fixed-seed parity test against pre-optimization output; full existing
robust-regression test suite.

### Phase 4: Re-benchmark and audit update
Re-run the Benchmark section's measurements; if `skip_brt="ci"` can be
lifted, update `path_audits_source.R:34` and regenerate `path_audits.html`.

## Acceptance Criteria

- Profiling breakdown is produced and documented before any optimization
  code is written.
- `compute_robust_rand_bootstrap_parallel_cpp` output matches
  pre-optimization output within `1e-10` on a fixed-seed test.
- Full existing robust-regression test suite (5 files) plus the new parity
  test all pass.
- Measured wall-clock improvement for `compute_fast_rand_bootstrap_distr` at
  realistic B is reported with actual numbers (no fixed target multiplier
  committed upfront).
- `estimate_only = TRUE` continues to skip all SE/vcov computation and
  allocation, verified by the profiling breakdown showing no regression in
  that path's allocation footprint.
