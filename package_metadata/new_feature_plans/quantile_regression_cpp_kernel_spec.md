# Quantile Regression C++ Kernel Spec

Generated: 2026-07-29

## Scope

This spec defines an implementation plan for a C++ kernel that fits quantile
regression (minimizes the pinball/check loss) with exact parity to
`quantreg::rq.fit(method = "br")` (point estimate) and
`quantreg:::summary.rq(se = "nid")` / `se = "iid"` (standard errors), so that
every class currently hard-dependent on `quantreg` gets a `use_rcpp` fast
path.

This is the largest of the four gaps found in the 2026-07-28 model-fitting
audit: unlike the ordinal-GEE, NegBin-weighted, and robust-regression gaps
(which each reused substantial existing C++/statistical infrastructure), no
quantile-regression code exists anywhere in `EDI/src/` today. The point
estimate requires a genuine linear-program (simplex) solver, not a
reformulation of existing gradient-based fitting code.

Related documents:

- [ordinal_gee_cpp_kernel_spec.md](ordinal_gee_cpp_kernel_spec.md) and
  [robust_regression_perf_optimization_spec.md](robust_regression_perf_optimization_spec.md) —
  companion specs from the same 2026-07-28 audit; establish the
  exact-parity-with-reference-implementation precedent this spec follows.
- [../package_tests/path_audits_source.R](../package_tests/path_audits_source.R)
  line 29 — `InferenceContinKKQuantileRegrOneLik` audit row: "bootstrap CI
  mean 109.0s / max 9434.3s at n=98 slow" — the concrete performance problem
  this spec addresses.

Affected classes (all confirmed to route through the same two `quantreg`
entry points — `rq.fit(x, y, tau, method = "br")` for `estimate_only = TRUE`,
`rq(formula, data, tau[, weights])` otherwise):

- `InferenceContinQuantileRegr` (`inference_continuous_quantile_regr.R`)
- `InferencePropQuantileRegr` (`inference_proportion_quantile_regr.R`,
  fits on a logit-transformed response — same LP, transform happens in R)
- `InferenceAbstractKKQuantileRegrOneLik` (`inference_all_KK_quantile_regr_one_lik_abstract.R`,
  parent of `InferenceContinKKQuantileRegrOneLik` /
  `InferencePropKKQuantileRegrOneLik`) — stacked matched-pair + reservoir
  design, same shape as the robust-regression combined-likelihood classes
- `InferenceAbstractKKQuantileRegrIVWC` (`inference_all_KK_quantile_regr_ivwc_abstract.R`) —
  separate matched-differences and reservoir fits combined via
  inverse-variance weighting, same shape as the robust-regression IVWC class

## Motivation

Every one of these classes calls into `quantreg` (an R/Fortran package) for
every single model fit — the point estimate, every asymptotic CI/p-value,
every randomization-inference replicate, and every bootstrap-weighted
refit. `path_audits_source.R:29` already documents this as a measured
performance problem (up to 9434s for a single bootstrap CI at n=98). This is
the last major inference family with a fully external, un-accelerated
dependency.

## Non-Goals

- Do not implement `quantreg`'s interior-point methods (`"fn"`, `"pfn"`,
  `"sfn"`) used for very large `n`. KK-design sample sizes are small to
  moderate (matched pairs + reservoir); `method = "br"` (Barrodale-Roberts
  simplex) is what this codebase already calls exclusively
  (`getFromNamespace("rq.fit", "quantreg")(..., method = "br")`), and is the
  correct target to replicate.
- Do not implement `se = "rank"` (the default for `n < 1001` in
  `quantreg:::summary.rq` when `se` is unspecified) — this codebase always
  passes `se` explicitly (`"nid"` then `"iid"` fallback,
  `.extract_se_from_rq_fit` in `other_helpers.R:624`), so `"rank"` is never
  reached and is out of scope.
- Do not implement `se = "boot"` or `se = "ker"` — not used by this
  codebase's `.extract_se_from_rq_fit`.
- Do not remove `quantreg` as a dependency or the `use_rcpp = FALSE`
  fallback path. `quantreg` remains the parity-test oracle and fallback,
  matching the pattern used for every other `use_rcpp`-gated class.
- Do not change the logit transform in `InferencePropQuantileRegr` or the
  matched-pair/reservoir design-stacking logic in the KK classes — those
  stay in R; only the underlying LP solve and SE computation move to C++.

## Statistical Background

Quantile regression at level `tau` minimizes the pinball (check) loss:

```text
minimize_beta  sum_i rho_tau(y_i - x_i' * beta)
rho_tau(u) = u * (tau - 1{u < 0})
```

This is equivalent to a linear program in the residual decomposition
`u_i = u_i^+ - u_i^-`, `u_i^+, u_i^- >= 0`:

```text
minimize            tau * sum(u^+) + (1 - tau) * sum(u^-)
subject to    X*beta + u^+ - u^- = y
              u^+, u^- >= 0
```

`quantreg::rq.fit(method = "br")` solves this via the Barrodale & Roberts
(1973) simplex algorithm as extended to regression quantiles by Koenker &
d'Orey (1987) — a specialized simplex method that exploits the piecewise-linear
structure of the pinball loss (pivoting through basic feasible solutions,
each of which exactly interpolates `p` observations, i.e. has `p` zero
residuals). This is implemented in `quantreg` as Fortran routine `rqbr.f`.

## Algorithm — Point Estimate

**Primary implementation reference: `quantreg`'s own `rqbr.f` Fortran
source**, not a from-memory description of the simplex method. Porting the
actual algorithm (rather than writing an independent simplex
implementation) is what guarantees exact numerical parity — this LP can have
non-unique or degenerate solutions at ties, and only matching the reference
implementation's specific pivoting/tie-breaking rules guarantees the C++
kernel picks the same vertex `quantreg` does.

Implementation steps:

1. Port `rqbr.f`'s simplex logic to C++ using `Eigen` for the linear algebra
   (consistent with every other kernel in `EDI/src/`), preserving its pivot
   selection and tie-breaking behavior exactly.
2. Validate against `quantreg::rq.fit(x, y, tau, method = "br")` directly
   (no R6 class involved yet) across many random `(X, y, tau)` triples,
   including degenerate/tied cases, before any R integration.

## Algorithm — Standard Errors

Both formulas below are transcribed directly from `quantreg:::summary.rq`
(confirmed against the installed package source, not recalled from memory)
and must be replicated exactly.

**`se = "iid"`:**

```r
xxinv = solve(t(x) %*% x)                     # via QR/backsolve of x
pz = sum(abs(resid) < eps)
h = max(p + 1, ceiling(n * bandwidth.rq(tau, n, hs = TRUE)))
ir = (pz + 1):(h + pz + 1)
ord.resid = sort(resid[order(abs(resid))][ir])
xt = ir / (n - p)
sparsity = rq(ord.resid ~ xt)$coef[2]          # one more small LP solve
cov = sparsity^2 * xxinv * tau * (1 - tau)
se = sqrt(diag(cov))
```

**`se = "nid"` (primary; the "Powell sandwich"):**

```r
h = bandwidth.rq(tau, n, hs = TRUE)            # Hall-Sheather bandwidth
while (tau - h < 0 || tau + h > 1) h = h / 2
bhi = rq.fit(x, y, tau = tau + h, method = "br")$coef   # 2 more LP solves
blo = rq.fit(x, y, tau = tau - h, method = "br")$coef
dyhat = x %*% (bhi - blo)
f = pmax(0, (2 * h) / (dyhat - eps))           # per-observation local sparsity
fxxinv = solve of (X' diag(f) X) via QR of sqrt(f)*X, then fxxinv %*% t(fxxinv)
cov = tau * (1 - tau) * fxxinv %*% (t(x) %*% x) %*% fxxinv
se = sqrt(diag(cov))
```

with Hall-Sheather bandwidth:

```r
bandwidth.rq(tau, n, hs = TRUE) =
  n^(-1/3) * qnorm(0.975)^(2/3) * ((1.5 * dnorm(qnorm(tau))^2) / (2 * qnorm(tau)^2 + 1))^(1/3)
```

Note `"nid"` requires **two additional LP solves** (at `tau + h` and
`tau - h`) beyond the point estimate — the C++ kernel's SE path is not a
closed-form add-on, it re-invokes the simplex solver twice more. This must
be reflected in the C++ interface (see below) and in perf expectations: the
`nid` path costs roughly 3x one LP solve, not 1x.

`.extract_se_from_rq_fit` (`other_helpers.R:624`) tries `"nid"` first and
falls back to `"iid"` if the result is non-finite, non-positive, or exceeds
`EDI_SEPARATION_THRESHOLD`. The C++ kernel must expose both and let the R
wrapper apply the same fallback logic, preserving current behavior exactly.

## C++ Interface

New file `EDI/src/fast_quantile_regression.cpp`:

```cpp
// [[Rcpp::export]]
List fast_rq_fit_br_cpp(SEXP X_r, SEXP y_r, double tau);
// Point estimate only — mirrors rq.fit(x, y, tau, method="br").
// Returns: list(coefficients, residuals, converged).

// [[Rcpp::export]]
List fast_rq_fit_br_weighted_cpp(SEXP X_r, SEXP y_r, SEXP weights_r, double tau);
// Weighted point estimate — mirrors rq(y ~ ., data, tau, weights)'s
// underlying weighted rq.fit call.

// [[Rcpp::export]]
List fast_rq_with_se_cpp(SEXP X_r, SEXP y_r, double tau, std::string se = "nid", bool hs = true);
// Point estimate + standard errors. Internally solves the LP up to 3 times
// for se="nid" (tau, tau+h, tau-h) or up to 2 times for se="iid" (tau, plus
// the small sparsity regression). Returns:
// list(coefficients, se, cov, sparsity_or_scale, se_method_used, converged)
```

`fast_rq_with_se_cpp` takes `se` as an explicit argument rather than always
computing both, so the R wrapper controls the nid-then-iid-fallback sequence
exactly as `.extract_se_from_rq_fit` does today: call once with `se="nid"`;
only call again with `se="iid"` if that result is non-finite, non-positive,
or exceeds `EDI_SEPARATION_THRESHOLD`. This keeps the common case (a usable
`nid` result) to a single call.

## R Integration

Each affected class gets `use_rcpp = TRUE` (default), following the
established pattern:

- `InferenceContinQuantileRegr` / `InferencePropQuantileRegr`: `fit_quantile_model`'s
  `estimate_only = TRUE` branch calls `fast_rq_fit_br_cpp`; the `FALSE`
  branch calls `fast_rq_with_se_cpp` with `se = "nid"`, falling back to
  `se = "iid"` on a bad result exactly as `.extract_se_from_rq_fit` does.
  `compute_estimate_with_bootstrap_weights` calls
  `fast_rq_fit_br_weighted_cpp`.
- `InferenceAbstractKKQuantileRegrOneLik`: same split, applied to the
  stacked matched-pair-difference + reservoir design (mirrors
  `InferenceContinKKRobustRegrOneLik`'s `fit_combined`/`fit_weighted_combined`
  structure from the robust-regression wiring).
- `InferenceAbstractKKQuantileRegrIVWC`: same split, applied separately to
  the matched-differences-only fit (`yd ~ 1`) and the reservoir fit
  (`yr ~ .`), combined via the existing inverse-variance weighting logic
  (mirrors `InferenceContinKKRobustRegrIVWC`'s `robust_for_matched_pairs`/
  `robust_for_reservoir` structure).
- `!use_rcpp` paths are untouched — `quantreg` stays exactly as it is today
  as the fallback and parity-test oracle.

## Testing

- Direct-solver parity test (no R6 involved): `fast_rq_fit_br_cpp` vs
  `quantreg::rq.fit(x, y, tau, method="br")` across many random `(X, y,
  tau)` triples, including near-degenerate/tied-residual cases, at tight
  tolerance (this is an exact LP, not an iterative approximation — mismatches
  beyond floating-point noise indicate a real bug, not an expected
  algorithmic difference).
- `fast_rq_with_se_cpp(se="nid")` and `se="iid"` vs
  `quantreg:::summary.rq(fit, se="nid"/"iid")$coefficients[, "Std. Error"]`,
  same tight-tolerance standard.
- Per-class `use_rcpp = TRUE` vs `FALSE` parity tests
  (`compute_estimate`, `compute_asymp_confidence_interval`,
  `compute_asymp_two_sided_pval`) for all four affected classes, following
  the `compare_kk_gee_wrapper_paths` pattern from `test-kk-gee-parity.R`.
- Weighted-fit parity test for `compute_estimate_with_bootstrap_weights`
  across all four classes.

## Perf / SIMD

Per the precedent set in `robust_regression_perf_optimization_spec.md`:
**profile before committing to SIMD work.** The Barrodale-Roberts simplex is
inherently sequential (pivot-by-pivot), so the primary perf lever is
unlikely to be intra-solve vectorization — it is more likely to be:

1. Avoiding the ~3x LP-solve cost of `se="nid"` where it isn't needed
   (`estimate_only = TRUE` paths, already handled by the interface split
   above).
2. The same per-call allocation-overhead concern flagged in the
   robust-regression spec, for classes that invoke this kernel repeatedly
   inside randomization-inference or bootstrap loops.
3. Any genuinely hot, vectorizable inner loop the profile identifies inside
   a single simplex pivot (e.g. ratio-test scans over residuals) — only
   pursue explicit SIMD here if profiling justifies it, consistent with the
   non-goal framing in the robust-regression spec.

## Benchmark

Before/after wall-clock for `path_audits_source.R:29`'s documented slow
path (`InferenceContinKKQuantileRegrOneLik`'s bootstrap CI at n≈98) plus a
smaller and larger n for context. Report actual numbers; no fixed target
multiplier committed upfront.

## Risks / Open Questions

- **Degenerate/tied-residual LP solutions.** The pinball-loss LP can have a
  non-unique optimal vertex at ties; exact parity depends on matching
  `rqbr.f`'s specific tie-breaking, not just reaching *an* optimal
  solution.
- **`nid`'s extra LP solves compound any LP-solver bugs.** Since `nid`
  re-invokes the point-estimate solver twice more (at `tau ± h`), any
  correctness issue in the core solver shows up three times as often in the
  SE path as in the point-estimate path — the direct-solver parity test
  above should be exhaustive before building the SE layer on top of it.
- **Weighted LP formulation.** Need to confirm during implementation exactly
  how `quantreg::rq(weights=...)` transforms the LP (likely `tau/(1-tau)`
  weighted per-observation costs analogous to the weighted-least-squares
  case, but must be verified against `quantreg` source rather than assumed).

## Implementation Phases

### Phase 1: Core LP solver
- [ ] TODO-1: Port `rqbr.f` to `fast_rq_fit_br_cpp`. Validate against
`quantreg::rq.fit(method="br")` directly, including degenerate cases.

### Phase 2: Weighted variant
- [ ] TODO-2: Implement `fast_rq_fit_br_weighted_cpp`, validated against weighted
`quantreg::rq`.

### Phase 3: Standard errors
- [ ] TODO-3: Implement `fast_rq_with_se_cpp` (`nid` and `iid`), validated against
`quantreg:::summary.rq`.

### Phase 4: R integration
- [ ] TODO-4: Wire `use_rcpp` through all four affected classes (`InferenceContinQuantileRegr`,
`InferencePropQuantileRegr`, `InferenceAbstractKKQuantileRegrOneLik`,
`InferenceAbstractKKQuantileRegrIVWC`).

### Phase 5: Test suite, profiling, and perf pass
- [ ] TODO-5: Full parity test suite per "Testing". Profile per "Perf / SIMD" before
implementing any perf-specific changes.

### Phase 6: Documentation and audit update
- [ ] TODO-6: Update `path_audits_source.R:29` (and the other three affected rows) to
record `use_rcpp` support and the measured speedup; regenerate
`path_audits.html`.

## Acceptance Criteria

- `fast_rq_fit_br_cpp` / `_weighted_cpp` match `quantreg::rq.fit`/`rq`
  point estimates within floating-point tolerance across the full test
  matrix in "Testing", including degenerate cases.
- `fast_rq_with_se_cpp` matches `quantreg:::summary.rq`'s `nid` and `iid`
  standard errors within floating-point tolerance.
- All four affected classes produce matching `compute_estimate`,
  `compute_asymp_confidence_interval`, `compute_asymp_two_sided_pval`, and
  bootstrap-weighted-estimate results under `use_rcpp = TRUE` vs `FALSE`.
- `quantreg` remains a `Suggests` dependency; no existing
  `use_rcpp = FALSE` behavior changes.
- Measured wall-clock improvement is reported for `path_audits_source.R:29`'s
  documented slow case; if large enough, that row's slow-path flag is
  reconsidered as part of Phase 6.
