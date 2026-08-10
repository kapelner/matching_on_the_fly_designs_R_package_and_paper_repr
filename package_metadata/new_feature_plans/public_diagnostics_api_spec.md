# Public Diagnostics API Implementation Spec

Generated: 2026-07-27

## Scope

This spec defines a public, opt-in diagnostics API for the primary inference
calls:

```r
compute_estimate()
compute_pval()
compute_ci()
```

The new API adds diagnostic variants:

```r
compute_estimate_details()
compute_pval_details()
compute_ci_details()
```

The existing scalar API must remain the fast path. It continues to return the
current value shape: a numeric estimate, p-value, confidence interval, or
`NA`/`c(NA_real_, NA_real_)` when the target is not estimable.

The debug API returns a structured object describing what happened during the
same logical computation: value, status, typed non-estimability reason,
timings, selected method, and collected diagnostics.

This feature should be implemented together with, or shortly after,
[optimizer_diagnostics_report.md](optimizer_diagnostics_report.md). That report
defines the lower-level solver/fit diagnostics layer. This spec defines the
public R-facing envelope that exposes those diagnostics without slowing normal
users.

## Non-Goals

- Do not change return values of `compute_estimate()`, `compute_pval()`, or
  `compute_ci()`.
- Do not make the fast path compute expensive diagnostics.
- Do not promise "all internals." The public debug object is a stable,
  curated diagnostic schema.
- Do not replace existing lower-level `debug = TRUE` options for bootstrap or
  randomization distribution helpers in the first pass. Those can feed into
  this API later.

## Naming

Use the exact public method names requested for now:

```r
compute_estimate_details()
compute_pval_details()
compute_ci_details()
```

Long term, other aliases (e.g. `diagnose_estimate()`) may be added, but the
initial implementation should keep one official name set to avoid documenting
multiple contracts.

## Core Design Rule

Shared implementation must share orchestration, not instrumentation.

The fast path must pass `diagnostics = FALSE` all the way down. Lower layers
must not compute extra QR decompositions, eigendecompositions, refits,
per-replicate traces, or root-search traces unless diagnostics are explicitly
requested.

Recommended shape:

```r
compute_estimate = function(...) {
  private$compute_estimate_core(..., diagnostics = FALSE)$value
}

compute_estimate_details = function(...) {
  private$compute_estimate_core(..., diagnostics = TRUE)
}
```

The core function may return the full debug object internally in both cases,
but when `diagnostics = FALSE` the object must contain only fields already
available at no extra cost. Avoid building large nested structures on the fast
path if allocation overhead is measurable.

## Diagnostic Object Contract

Return an S3 object with class:

```r
c("EDIInferenceDebugResult", "list")
```

Required top-level fields:

```r
list(
  value = <numeric result>,
  status = "ok" | "nonestimable" | "error",
  target = "estimate" | "pval" | "ci",
  method = <character scalar>,
  inference_class = <character scalar>,
  response_type = <character scalar>,
  design_type = <character scalar>,
  stage = <character scalar or NULL>,
  reason = <character scalar or NULL>,
  message = <character scalar or NULL>,
  timings = list(total_sec = <numeric>),
  diagnostics = list(),
  warnings = character(),
  errors = character()
)
```

`value` must match the public scalar API for the same object and arguments.
For non-estimable cases:

- estimate: `NA_real_`
- p-value: `NA_real_`
- CI: `c(NA_real_, NA_real_)` with normal CI names where applicable

`status` meanings:

- `ok`: the returned value is finite and usable.
- `nonestimable`: the path completed and returned a typed non-estimable state.
- `error`: an unexpected exception occurred. The debug API may catch and return
  this as a structured debug result, but the existing scalar API behavior
  should remain unchanged unless explicitly changed elsewhere.

`stage` should use the package's existing stages where possible:

```r
"estimate"
"se"
"pval"
"ci"
"bootstrap"
"bayesian_bootstrap"
"parametric_bootstrap"
"jackknife"
"randomization"
"ci_inversion"
"likelihood_null_fit"
"likelihood_full_fit"
"solver"
```

`reason` must be a stable typed string, not a free-form sentence. Existing
`cache_nonestimable_estimate()` and `cache_nonestimable_se()` reasons should
be reused where possible.

## Method-Specific Required Fields

### `compute_estimate_details()`

`estimate_only` is an explicit input parameter here, mirroring the existing
public `compute_estimate(estimate_only = FALSE)` signature -- the sole
argument on every one of the ~85 `compute_estimate()` overrides across the
package, with no exceptions -- rather than something only discovered after
the fact:

```r
compute_estimate_details = function(estimate_only = FALSE, diagnostic_level = c("standard", "extended", "trace"))
```

The `diagnostic_level` argument is the only addition beyond what
`compute_estimate()` already takes; see Lazy Diagnostics Levels below.

Required `diagnostics` subfields when available without extra work:

```r
diagnostics = list(
  estimate_only = TRUE | FALSE,
  converged = TRUE | FALSE | NA,
  iterations = <integer or NA>,
  coefficient_abs_max = <numeric or NA>,
  standard_error = <numeric or NA>,
  nonestimable_state = <list or NULL>,
  fit = <optimizer diagnostics or NULL>
)
```

`diagnostics$estimate_only` echoes the effective value actually used for the
underlying fit this snapshot describes -- normally the caller's own
`estimate_only` argument above, passed straight through the same
`private$shared(estimate_only = ...)` path `compute_estimate()` already uses.
It can still differ from the caller's input when, at diagnostic level
`"standard"`, an already-cached fit is reused instead of a fresh one (see
Lazy Diagnostics Levels) -- in that case it reports the mode the reused fit
actually ran in, so the caller can tell why `standard_error` or `fit` came
back `NA` even though they passed `estimate_only = FALSE`.

Setting `estimate_only = TRUE` still computes and returns `value` -- it only
concerns which of the diagnostics subfields above are cheaply available;
`converged`, `iterations`, and `coefficient_abs_max` are populated by any fit
mode, while `standard_error` and `fit$information`-style fields require
`estimate_only = FALSE`.

Do not compute rank, condition number, or eigenspectrum unless diagnostics mode
explicitly asks for them and the lower layer advertises support.

### `compute_pval_details()`

Required arguments mirrored in the result:

```r
parameters = list(
  delta = <numeric>,
  testing_type = <character>,
  alternative = "two_sided"
)
```

Required `diagnostics` subfields where applicable:

```r
diagnostics = list(
  estimate = <numeric or NA>,
  standard_error = <numeric or NA>,
  test_statistic = <numeric or NA>,
  reference_distribution = <character or NULL>,
  full_fit = <fit diagnostics or NULL>,
  null_fit = <fit diagnostics or NULL>,
  bootstrap = <bootstrap diagnostics or NULL>,
  randomization = <randomization diagnostics or NULL>
)
```

Bootstrap-style p-value diagnostics should include:

```r
list(
  n_requested = <integer>,
  n_attempted = <integer>,
  n_success = <integer>,
  n_finite = <integer>,
  n_nonfinite = <integer>,
  n_extreme = <integer>,
  min_required = <integer>,
  finite_fraction = <numeric>,
  reused_workers = TRUE | FALSE | NA,
  seed = <integer or NULL>
)
```

### `compute_ci_details()`

Required arguments mirrored in the result:

```r
parameters = list(
  alpha = <numeric>,
  testing_type = <character or NULL>,
  confidence_level = 1 - alpha
)
```

Required `diagnostics` subfields where applicable:

```r
diagnostics = list(
  estimate = <numeric or NA>,
  standard_error = <numeric or NA>,
  interval_method = <character>,
  fallback_used = TRUE | FALSE | NA,
  fallback_method = <character or NULL>,
  inversion = <list or NULL>,
  bootstrap = <list or NULL>
)
```

CI inversion diagnostics should include:

```r
list(
  root_engine = <character>,
  lower_seed = <numeric or NA>,
  upper_seed = <numeric or NA>,
  n_evaluations = <integer or NA>,
  bracket_found = TRUE | FALSE | NA,
  reject_at_estimate = TRUE | FALSE | NA,
  p_at_estimate = <numeric or NA>,
  validation_failure = <character or NULL>
)
```

## Lazy Diagnostics Levels

Use a `diagnostic_level` argument on debug methods:

```r
compute_estimate_details(..., diagnostic_level = c("standard", "extended", "trace"))
compute_pval_details(..., diagnostic_level = c("standard", "extended", "trace"))
compute_ci_details(..., diagnostic_level = c("standard", "extended", "trace"))
```

Default: `"standard"`.

Levels:

- `standard`: return values already computed by the normal path plus typed
  non-estimability state and timings. No extra model fits or expensive matrix
  decompositions.
- `extended`: allow one-off extra diagnostics such as QR rank checks or
  condition-number checks when the class advertises support and the cost is
  documented.
- `trace`: allow large per-replicate/per-root-evaluation traces. This is for
  development and should not be used by default in simulation sweeps.

The scalar fast API must never call `diagnostic_level = "extended"` or
`"trace"` internally.

### Why A Tiered Level Instead Of A Per-Diagnostic Flags List

An alternative design would let the caller pass a list enumerating exactly
which diagnostics to compute, e.g.
`diagnostics = list(condition_number = TRUE, quadrature_adequacy = FALSE, qr_rank = TRUE, ...)`.
This was rejected for the public contract:

- The set of *possible* diagnostics differs by family, and most flags would
  not apply to most families. `condition_number` only means something for
  likelihood-backed asymptotic families with an information matrix;
  `quadrature_adequacy` only exists for GH-quadrature GLMM paths; bootstrap
  finite-count diagnostics do not apply to closed-form Wald paths at all. A
  flags list forces every family to either silently ignore flags it does not
  support (confusing -- you asked for something and got nothing, with no
  signal why) or error on them (brittle -- breaks exactly the kind of generic
  sweep this package already relies on, e.g. `path_audits.html` and
  comprehensive tests calling the same method uniformly across dozens of
  heterogeneous classes).
- A named-flags API is functionally "all internals, gated behind booleans,"
  which contradicts the Non-Goals section above ("Do not promise 'all
  internals.' The public debug object is a stable, curated diagnostic
  schema"). A cost tier lets each family decide locally what `"standard"` vs
  `"extended"` vs `"trace"` means for itself without the public contract ever
  having to enumerate a global diagnostic vocabulary -- one that would only
  grow as [optimizer_diagnostics_report.md](optimizer_diagnostics_report.md)'s
  taxonomy (already 8 distinct failure categories) expands.

This is a real, acknowledged tradeoff, not a free win: the tiered design gives
up fine-grained cost control. A caller who wants only the condition number
and nothing else cannot ask for that in isolation -- they get whatever else
`"extended"` happens to bundle for that family, and pay for all of it. If a
specific diagnostic later turns out to be expensive enough that bundling it
into `"extended"` becomes a real cost problem for callers who don't want it,
that is the point to revisit this decision (e.g. an optional, family-specific
override layered on top of the tier), not to redesign the whole contract
around per-flag toggles up front.

## Integration With Optimizer Diagnostics

[optimizer_diagnostics_report.md](optimizer_diagnostics_report.md) proposes
native and R-level diagnostics for:

- iteration cap
- terminal gradient norm
- near-singular information
- separation or boundary divergence
- huge finite standard errors from an ill-conditioned information matrix
- GLMM variance-component boundary
- quadrature adequacy
- derivative-free optimizer convergence codes
- extreme bootstrap or LR statistics

The public details API should preserve the distinction between these cases.
Recent comprehensive-test failures in `InferenceOrdinalStereotypeLogitRegr`,
`InferenceOrdinalKKCondAdjCatLogitRegr`, and
`InferenceIncidKKCondLogitOneLik` showed two different shapes: extreme finite
point estimates, and moderate finite point estimates paired with enormous
finite SEs. Both should produce typed diagnostics; they should not be collapsed
into a generic "separation" label unless the coefficient path itself is the
evidence.

When that layer lands, debug methods should embed its output under:

```r
diagnostics$fit
diagnostics$full_fit
diagnostics$null_fit
diagnostics$replicate_fits
```

Expected fit diagnostic schema:

```r
list(
  converged = TRUE | FALSE | NA,
  hit_iteration_cap = TRUE | FALSE | NA,
  iterations = <integer or NA>,
  gradient_norm = <numeric or NA>,
  information_min_eigenvalue = <numeric or NA>,
  information_condition_number = <numeric or NA>,
  coefficient_abs_max = <numeric or NA>,
  boundary_hit = TRUE | FALSE | NA,
  boundary_type = <character or NULL>,
  optimizer = <character or NULL>,
  optimizer_code = <integer or NA>,
  diagnostic_category = <character or NULL>,
  diagnostic_severity = "ok" | "warning" | "failure" | NA
)
```

If optimizer diagnostics are not yet available for a family, set the relevant
field to `NULL` or `NA`, not by doing expensive substitute computation.

## Internal Implementation Plan

### Phase 1: Public wrapper and result object

- [ ] TODO-1: Add an `EDIInferenceDebugResult` constructor/helper.
- [ ] TODO-2: Add print and summary methods:
   - `print.EDIInferenceDebugResult`
   - `as.data.frame.EDIInferenceDebugResult`
- [ ] TODO-3: Add base public methods in the highest shared inference class that can
   safely expose them.
- [ ] TODO-4: The default implementations should call existing public methods and wrap
   their results with minimal diagnostics. This makes the API available
   package-wide without changing class internals.

Phase 1 default behavior is intentionally shallow but safe.

### Phase 2: Core-path integration

Move high-traffic families from wrapper-only diagnostics to shared core
functions:

```r
private$compute_estimate_core(diagnostics = FALSE, diagnostic_level = "standard")
private$compute_pval_core(diagnostics = FALSE, diagnostic_level = "standard")
private$compute_ci_core(diagnostics = FALSE, diagnostic_level = "standard")
```

Do this family-by-family. Do not attempt a single large rewrite across all
classes.

Priority targets:

- [ ] TODO-5: Move likelihood-backed `InferenceAsympLik` classes to shared core-path diagnostics.
- [ ] TODO-6: Move `InferenceParamBootstrap` and bootstrap-calibrated LR paths to shared core-path diagnostics.
- [ ] TODO-7: Move nonparametric bootstrap and Bayesian bootstrap methods to shared core-path diagnostics.
- [ ] TODO-8: Move randomization and jackknife paths to shared core-path diagnostics.

### Phase 3: Optimizer diagnostics integration

After the optimizer diagnostics layer is implemented:

- [ ] TODO-9: Thread native diagnostics into R fit objects.
- [ ] TODO-10: Store last-fit diagnostics in a consistent private cache.
- [ ] TODO-11: Expose `get_last_fit_diagnostics()` as planned by
   [optimizer_diagnostics_report.md](optimizer_diagnostics_report.md).
- [ ] TODO-12: Have `compute_*_details()` pull from those caches without recomputing.

### Phase 4: Audit/report integration

Use debug results to improve `path_audits.html` and comprehensive tests:

- [ ] TODO-13: Distinguish numeric success, typed non-estimable, and unexpected error.
- [ ] TODO-14: Aggregate `reason` by class/method.
- [ ] TODO-15: Identify dominant failure mechanisms for 0%, <1%, <5%, and <25% cells.
- [ ] TODO-16: Produce a low-estimability hardening report directly from debug output.

This should complement [path_audit_hardening_report.md](../package_tests/path_audit_hardening_report.md).

## Fast-Path Performance Requirements

The following are hard requirements:

1. `compute_estimate()`, `compute_pval()`, and `compute_ci()` must not allocate
   large diagnostic structures.
2. They must not compute QR, eigenvalues, condition numbers, refits, or
   replicate traces solely for diagnostics.
3. They must not call debug wrappers internally.
4. Any shared core must branch before expensive instrumentation.
5. Benchmarks for representative fast paths must show no meaningful regression.

Suggested benchmark targets:

- closed-form/simple mean difference
- OLS
- logistic/probit GLM
- Poisson/count GLM
- one GLMM path
- one bootstrap path

Acceptance threshold: median runtime change on scalar calls should be within
measurement noise, with no systematic allocation increase visible in `Rprofmem`
or comparable tooling.

## Error Handling

Debug methods should catch unexpected errors and return:

```r
status = "error"
value = NA_real_ # or c(NA_real_, NA_real_) for CI
stage = <best known stage>
reason = "unexpected_error"
message = conditionMessage(e)
errors = conditionMessage(e)
```

This is a diagnostic API difference from the scalar API. Scalar methods should
keep their current behavior unless a separate API decision changes it.

## Public Documentation

Document that:

- scalar methods are for routine analysis and simulation sweeps
- debug methods are for investigating `NA`, non-estimable output, numerical
  fragility, and audit failures
- debug methods may be slower and may return larger objects
- not every diagnostic field is available for every class
- missing diagnostic fields mean "not collected" or "not applicable", not
  necessarily "no problem"

## Minimal Example

```r
res = inf$compute_pval_details(delta = 0)

res$value
res$status
res$stage
res$reason
res$diagnostics$fit$converged
res$diagnostics$bootstrap$n_finite
```

Example non-estimable output:

```r
list(
  value = NA_real_,
  status = "nonestimable",
  target = "pval",
  method = "lik_ratio_bartlett_approx",
  stage = "se",
  reason = "lik_ratio_bartlett_approx_test_unavailable",
  diagnostics = list(
    full_fit = list(converged = TRUE, iterations = 8),
    null_fit = list(converged = FALSE, diagnostic_category = "separation"),
    bootstrap = list(n_requested = 151, n_finite = 17, min_required = 31)
  )
)
```

## Acceptance Criteria

- `compute_estimate_details()`, `compute_pval_details()`, and
  `compute_ci_details()` exist for all public inference objects.
- Existing scalar methods return exactly the same values as before.
- Debug `value` matches the scalar return for the same arguments.
- Debug objects carry stable `status`, `stage`, and `reason` fields.
- At least one likelihood-backed family includes optimizer diagnostics once the
  optimizer diagnostics layer lands.
- At least one bootstrap path includes finite-count diagnostics.
- Fast-path benchmark shows no meaningful regression.
