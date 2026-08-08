# PRW Subsampling Implementation Spec

Generated: 2026-07-27

## Scope

This spec defines an implementation plan for Politis/Romano/Wolf-style
subsampling inference in EDI.

PRW subsampling is a resampling method based on drawing subsamples of size
`b < n` without replacement, recomputing the statistic on each subsample, and
using the empirical distribution of the properly centered and scaled subsample
statistics to approximate the sampling distribution of the full-sample
estimator.

This is not the ordinary nonparametric bootstrap and should not be implemented
as a `bootstrap_type` option. The current bootstrap infrastructure draws
exchangeable units with replacement and usually draws `n` units. PRW
subsampling draws `b` units without replacement and uses different calibration.

Related implementation documents:

- [m_out_of_n_bootstrap_implementation_spec.md](m_out_of_n_bootstrap_implementation_spec.md)
- [public_diagnostics_api_spec.md](public_diagnostics_api_spec.md)
- [optimizer_diagnostics_report.md](optimizer_diagnostics_report.md)
- [../package_tests/path_audit_hardening_report.md](../package_tests/path_audit_hardening_report.md)

## Motivation

Subsampling gives EDI a principled fallback for paths where ordinary bootstrap
inference is questionable because of nonregularity, boundary behavior,
nonsmooth estimators, matching effects, or bootstrap refit instability.

It should be treated as a new inference family with explicit diagnostics, not
as a silent repair for failed bootstrap draws.

## Non-Goals

- Do not change the existing `compute_bootstrap_*()` behavior.
- Do not silently route ordinary bootstrap failures to subsampling.
- Do not implement m-out-of-n bootstrap in this spec. It is related but uses
  with-replacement samples and needs a separate spec.
- Do not claim universal validity. Subsampling still requires asymptotic
  assumptions and a defensible choice of subsample size.
- Do not make the scalar fast paths compute subsampling diagnostics unless a
  subsampling method is explicitly called.

## Public API

Add explicit public methods on inference objects that inherit the shared
resampling stack:

```r
approximate_subsampling_distribution_beta_hat_T(
  B = 501,
  b = NULL,
  show_progress = TRUE,
  debug = FALSE,
  subsampling_type = NULL,
  scaling = "sqrt_n",
  center = "full_estimate"
)

compute_subsampling_two_sided_pval(
  delta = 0,
  B = 501,
  b = NULL,
  type = "centered",
  show_progress = TRUE,
  min_number_usable_samples = 5L,
  subsampling_type = NULL,
  scaling = "sqrt_n"
)

compute_subsampling_confidence_interval(
  alpha = 0.05,
  B = 501,
  b = NULL,
  type = "basic",
  show_progress = TRUE,
  min_number_usable_samples = 5L,
  subsampling_type = NULL,
  scaling = "sqrt_n"
)
```

Do not overload `compute_bootstrap_two_sided_pval()` or
`compute_bootstrap_confidence_interval()` with a subsampling mode. The method
has different theory and users should see that difference at the API boundary.

## Theory Contract

For a full-sample estimate `theta_hat_n`, a subsample estimate
`theta_hat_b_star`, and a scaling sequence `tau_n`, PRW subsampling estimates
the distribution of:

```text
tau_n * (theta_hat_n - theta)
```

with the empirical distribution of:

```text
tau_b * (theta_hat_b_star - theta_hat_n)
```

where `b -> infinity` and `b / n -> 0` asymptotically.

For root-n regular estimators, the default scaling is:

```text
tau_n = sqrt(n)
tau_b = sqrt(b)
```

The confidence interval should use subsampling quantiles of the centered
subsample statistic:

```text
q_lo = quantile(tau_b * (theta_hat_b_star - theta_hat_n), alpha / 2)
q_hi = quantile(tau_b * (theta_hat_b_star - theta_hat_n), 1 - alpha / 2)

CI = c(theta_hat_n - q_hi / tau_n,
       theta_hat_n - q_lo / tau_n)
```

The centered two-sided p-value for `H0: theta = delta` should compare:

```text
T_n = tau_n * (theta_hat_n - delta)
```

to the centered subsampling distribution:

```text
T_b_star = tau_b * (theta_hat_b_star - theta_hat_n)
```

The default p-value should be:

```text
2 * min(mean(T_b_star <= T_n), mean(T_b_star >= T_n))
```

clamped to `[0, 1]`.

## Scaling Options

Initial implementation should support:

```r
scaling = "sqrt_n"
```

Internally define:

```r
subsampling_scaling_factor = function(size, scaling) {
  if (identical(scaling, "sqrt_n")) sqrt(size)
}
```

Later extensions may add:

```r
scaling = list(rate_exponent = 0.5)
scaling = function(size) ...
```

but the first implementation should keep this narrow. Incorrect scaling can
invalidate inference, so the default should be explicit and documented.

## Choosing `b`

The first implementation should accept either:

```r
b = <integer>
b = NULL
b = list(method = "minimum_volatility", ...)
```

The implementation should support two levels of automatic selection:

1. a cheap deterministic rule for simple calls
2. a data-driven minimum-volatility selector for production use

The deterministic fallback rule is:

```r
b = floor(n_units ^ 0.7)
```

with guards:

```r
b >= max(5L, p_effective + 2L)
b <= floor(n_units / 2)
```

where `n_units` is the number of exchangeable subsampling units, not
necessarily the number of rows.

If the auto rule violates these bounds, return typed non-estimable output:

```text
reason = "subsampling_invalid_subsample_size"
stage  = "subsampling"
```

The debug API should report the chosen `b`, `n_units`, the rule used, and any
bound adjustments.

### Minimum-Volatility Selection

EDI should adapt the minimum-volatility approach already implemented in the PTE
package at `../PTE/PTE/R/select_optimal_m_prop.R`.

That PTE selector evaluates a grid of candidate exponent values, computes an
inference answer at each grid point, measures the local volatility of the
answer over a centered neighborhood, and chooses the interior grid point whose
answer is most stable. The same idea applies to PRW subsampling, replacing the
m-out-of-n bootstrap resample size `m` with the subsampling size `b`.

The selector logic should be shared with the m-out-of-n bootstrap implementation
described in
[m_out_of_n_bootstrap_implementation_spec.md](m_out_of_n_bootstrap_implementation_spec.md).
Both methods need to choose a tuning size that grows but remains small relative
to the full sample:

```text
m-out-of-n bootstrap: m -> infinity, m / n -> 0
PRW subsampling:      b -> infinity, b / n -> 0
```

The common selector should not know how to draw resamples. It should only
evaluate a grid and choose the least volatile eligible interior point:

```r
select_optimal_resample_size(
  size_grid,
  evaluator,
  objective = "ci_width",
  volatility_window = 3L,
  min_finite_fraction = 0.8
)
```

The evaluator is method-specific:

```r
evaluate_m_out_of_n_bootstrap_size(m, ...)
evaluate_subsampling_size(b, ...)
```

So the reusable part is minimum-volatility selection, not the resampling
implementation. m-out-of-n bootstrap samples with replacement; PRW subsampling
samples without replacement and uses its own operation contract.

Add a public calibration helper:

```r
select_optimal_b_subsampling(
  B = 251,
  alpha = 0.05,
  b_pow_of_n_grid = seq(0.5, 0.9, by = 0.05),
  b_grid = NULL,
  objective = "ci_width",
  target = "ci",
  volatility_window = 3L,
  subsampling_type = NULL,
  scaling = "sqrt_n",
  show_progress = TRUE
)
```

If `b_grid` is supplied, use it directly. Otherwise compute:

```r
b_grid = unique(floor(n_units ^ b_pow_of_n_grid))
```

and then apply the normal `b` validity guards. Grid points that fail validity
checks should be retained in the returned table with a typed failure reason,
but they are not eligible for selection.

The default exponent grid should avoid values too close to `1` because PRW
subsampling requires `b / n -> 0`. Values near `1` can be allowed explicitly as
diagnostic reference points, but should not be eligible by default.

Recommended first-pass defaults:

```r
b_pow_of_n_grid = seq(0.5, 0.9, by = 0.05)
volatility_window = 3L
B = 251
objective = "ci_width"
```

The first implementation should support only:

```r
objective = "ci_width"
target = "ci"
```

Additional objectives can be added after the CI-width selector is stable:

```r
objective = "pval_stability"
objective = "estimate_stability"
objective = "finite_fraction_penalized_ci_width"
```

For any objective, the selector should prefer stable interior regions over
edge optima. The objective value used for volatility must be recorded in the
returned `grid_table`.

The selector should compute, for each eligible `b`:

```r
estimate
ci_lower
ci_upper
ci_width
pval
objective_value
n_finite
finite_fraction
dominant_failure_reason
elapsed_sec
```

For the default CI-width objective, define local volatility as:

```r
sd(ci_width[(k - h):(k + h)], na.rm = FALSE)
```

where `h = (volatility_window - 1) / 2`, and only complete centered windows are
eligible. Boundary grid points should have `volatility = NA` and should not be
selectable, matching the PTE implementation. This avoids falsely selecting an
edge point only because it has a one-sided neighborhood.

The selected value is:

```r
k_star = eligible[which.min(volatility[eligible])]
b_optimal = b_grid[k_star]
```

Return an object:

```r
list(
  b_optimal = <integer>,
  b_pow_of_n_optimal = <numeric or NA>,
  grid_table = <data.frame>,
  objective = <character>,
  target = <character>,
  B = <integer>,
  alpha = <numeric>,
  volatility_window = <integer>,
  subsampling_type = <character or NULL>,
  scaling = <character>
)
```

with class:

```r
"EDISubsamplingBSelection"
```

### How `b = NULL` Should Use Selection

For the first implementation, `b = NULL` should use the deterministic
`floor(n_units ^ 0.7)` rule so that ordinary calls remain predictable and not
silently expensive.

After `select_optimal_b_subsampling()` is implemented and tested, support:

```r
b = list(method = "minimum_volatility")
```

This should run the selector with modest calibration settings, then run the
final requested inference at `b_optimal`. The final inference must be a fresh
run; the low-`B` calibration grid should not be treated as production
inference.

The debug API should record both stages:

```r
diagnostics$subsampling$b_selection
diagnostics$subsampling$final_run
```

### Sensitivity Helper

Add a related sensitivity helper:

```r
compute_subsampling_sensitivity(
  b_grid = NULL,
  b_pow_of_n_grid = seq(0.5, 0.9, by = 0.05),
  ...
)
```

This should return the same grid table as the selector, but without choosing a
single `b`. It is useful for documentation, reports, and debugging paths where
minimum volatility is ambiguous.

### Ambiguous Or Failed Selection

If all eligible grid points fail, return typed non-estimable output:

```text
reason = "subsampling_b_selection_failed"
stage  = "subsampling"
```

If multiple grid points have nearly identical minimum volatility, choose the
smallest `b` among ties by default. This keeps the selected value more aligned
with the `b / n -> 0` requirement. The tie rule should be recorded in
diagnostics.

If the selected grid point has an unacceptable finite fraction, return typed
non-estimable output rather than selecting the least-bad unstable value:

```text
reason = "subsampling_b_selection_too_few_finite_estimates"
stage  = "subsampling"
```

The minimum acceptable finite fraction should be configurable, with a default
of:

```r
min_finite_fraction = 0.8
```

## Design-Aware Subsampling Units

Subsampling must operate on exchangeable units, not blindly on rows.

Use a design-aware unit resolver analogous to jackknife unit resolution:

```r
resolve_subsampling_unit(unit = "auto")
```

Recommended unit mapping:

| Design family | Default subsampling unit |
|---|---|
| ordinary non-blocking designs | observation |
| cluster designs | cluster |
| blocked cluster designs | cluster within stratum |
| blocking designs | observation within block/stratum |
| matched-pair designs | pair |
| KK matched/reservoir designs | matched set plus reservoir observation |
| custom designs | observation, with a diagnostic warning |

For blocking designs, the first pass should preserve observed stratum
composition approximately by sampling within strata:

```text
b_h = round(b * n_h / n_units)
```

with deterministic correction so `sum_h b_h = b`.

If a stratum has too few units to draw the requested `b_h`, reduce `b_h` only
if the total `b` contract can still be met; otherwise return typed
non-estimable output:

```text
reason = "subsampling_stratum_too_small"
stage  = "subsampling"
```

For sequential designs, preserve the original order of sampled units in the
subsample object after drawing the unit set without replacement. This keeps
arrival-order semantics stable for estimators that inspect design history.

## Internal Architecture

Add a new resampling operation contract:

```r
subsampling = list(
  operation = "subsampling",
  draw_type = "unit_subsample_without_replacement",
  loader = "load_subsampling_draw_into_worker",
  estimator = "compute_subsampling_worker_estimate",
  cache_name = "subsampling_distr_cache",
  cache_key_method = "subsampling_cache_key"
)
```

Do not reuse `"non_param_boot"` because loaders, cache keys, draw metadata, and
distribution scaling differ.

Add a new abstract mixin/class in the resampling stack:

```r
InferenceSubsampling
```

Recommended inheritance placement:

```text
InferenceNonParamBootstrap
  -> InferenceRandBootstrap
  -> InferenceRandBootstrapCI
  -> InferenceBayesianBootstrap
  -> InferenceJackknife
```

The cleanest insertion point is alongside `InferenceNonParamBootstrap`, not
inside it. If the current inheritance chain makes that disruptive, implement
the first pass as private methods inside `InferenceNonParamBootstrap` but keep
the operation name, cache, and public methods explicitly subsampling-specific.

## Draw Representation

Represent each draw as:

```r
list(
  i_b = <integer row indices>,
  m_vec_b = <integer or NULL>,
  unit_ids = <integer unit ids>,
  unit_type = <character>,
  b = <integer>,
  n_units = <integer>,
  stratum_ids = <integer or NULL>
)
```

The row index vector `i_b` must contain no duplicate exchangeable units. It may
contain multiple rows per unit when the unit is a cluster, block, pair, or
matched set.

The draw generator should be private:

```r
private$subsampling_sample_indices(b, subsampling_type = NULL)
```

The default ordinary-design draw is:

```r
sample.int(n, b, replace = FALSE)
```

For performance, add a C++ helper only after the R implementation is correct
and profiled. The initial method can use base R sampling because the dominant
cost is usually refitting.

## Distribution Function

`approximate_subsampling_distribution_beta_hat_T()` should return a numeric
vector of raw subsample estimates when `debug = FALSE`, matching the current
bootstrap distribution helper style.

When `debug = TRUE`, return:

```r
list(
  values = <numeric raw theta_hat_b estimates>,
  centered_scaled_values = <numeric T_b_star values>,
  full_estimate = <numeric>,
  b = <integer>,
  n_units = <integer>,
  scaling = <character or list>,
  unit_type = <character>,
  errors = <list>,
  warnings = <list>,
  num_errors = <integer vector>,
  num_warnings = <integer vector>,
  prop_iterations_with_errors = <numeric>,
  prop_iterations_with_warnings = <numeric>,
  prop_illegal_values = <numeric>,
  finite_fraction = <numeric>
)
```

The scalar CI and p-value methods should compute centered/scaled values
internally even when the raw distribution helper returns raw estimates.

## Failure Handling

Subsampling should follow the same typed non-estimability discipline as the
diagnostics API.

Required reasons:

```text
subsampling_original_estimate_unavailable
subsampling_invalid_subsample_size
subsampling_too_few_finite_estimates
subsampling_stratum_too_small
subsampling_unit_structure_unavailable
subsampling_extreme_confidence_interval
```

Do not silently condition on successful subsamples when the finite fraction is
low. The method may drop non-finite subsample estimates only if:

```r
finite_count >= min_number_usable_samples
```

and diagnostics record:

```r
n_requested
n_attempted
n_success
n_finite
n_nonfinite
finite_fraction
failure_reasons
```

For the public scalar methods, default behavior should be:

- return `NA_real_` for p-values when too few finite subsamples exist
- return `c(NA_real_, NA_real_)` for intervals when too few finite subsamples
  exist
- cache a typed non-estimable state for the harness and debug API

## Diagnostics Integration

`compute_pval_debug()` and `compute_ci_debug()` from
[public_diagnostics_api_spec.md](public_diagnostics_api_spec.md) should expose
subsampling diagnostics under:

```r
diagnostics$subsampling
```

Required fields:

```r
list(
  n_units = <integer>,
  b = <integer>,
  B = <integer>,
  scaling = <character>,
  center = "full_estimate",
  unit_type = <character>,
  n_attempted = <integer>,
  n_finite = <integer>,
  finite_fraction = <numeric>,
  failure_reasons = <table/list>,
  q_lo = <numeric or NA>,
  q_hi = <numeric or NA>
)
```

If [optimizer_diagnostics_report.md](optimizer_diagnostics_report.md) is
implemented before this feature, failed subsample refits should preserve the
per-fit optimizer diagnostics in debug mode. The scalar subsampling API should
not collect full optimizer traces by default.

## Cache Keys

Subsampling cache keys must include:

```text
operation
B
b
subsampling_type
scaling
center
seed/draw identity if applicable
```

Do not share caches with ordinary bootstrap distributions. The raw values may
look similar, but the operation and calibration are different.

Selection/calibration cache keys must be separate from final-inference cache
keys and include:

```text
b_grid or b_pow_of_n_grid
selection_B
alpha
objective
target
volatility_window
min_finite_fraction
subsampling_type
scaling
```

Do not reuse a low-`B` selector run as the final reported subsampling
distribution unless the user explicitly requested the same `B`, `b`, and
arguments for final inference.

## Tests

Add focused tests before broad path-audit integration.

### Draw Tests

- ordinary design draws exactly `b` unique observations
- blocking design draws preserve approximate stratum allocation
- cluster design draws no duplicate clusters
- matched design draws intact pairs/matched sets
- sequential design subsample rows are returned in original arrival order

### Distribution Tests

- raw distribution length is `B`
- debug mode returns raw and centered/scaled distributions
- cache keys distinguish `b`, `B`, and `scaling`
- non-finite estimates are counted and typed

### Inference Tests

Use simple mean difference as the first correctness target:

- subsampling CI has ordered finite bounds in a benign dataset
- p-value is in `[0, 1]`
- larger `B` with fixed seed is deterministic
- invalid `b` returns typed non-estimable output where hardening is enabled

Then add one representative model-fit path:

- continuous OLS or incidence logistic regression with a benign dataset
- a sparse dataset that intentionally produces some failed subsamples

### `b` Selection Tests

- explicit `b_grid` is respected
- `b_pow_of_n_grid` maps to unique integer `b` values
- boundary grid points are not eligible when the volatility window is centered
- the minimum-volatility selector chooses the least volatile eligible interior
  point
- tie handling chooses the smallest `b`
- failed grid points remain in `grid_table` with typed reasons
- selector output is deterministic under a fixed seed
- final inference with `b = list(method = "minimum_volatility")` reruns at
  `b_optimal` instead of reusing the low-`B` selector distribution

### Harness Tests

Update `package_tests/path_audits_source.R` only after the method is implemented
and targeted tests pass. Add separate columns or a separate section for:

```text
subsample.pval
subsample.ci
```

Do not fold subsampling into existing bootstrap columns.

## Documentation

Document this as "subsampling" rather than "subsampling bootstrap" in the main
API names. In prose, mention Politis/Romano/Wolf and explain that it is related
to, but distinct from, the bootstrap.

Required documentation caveats:

- validity depends on `b -> infinity` and `b / n -> 0`
- the default `b = NULL` uses `floor(n_units ^ 0.7)`, a deterministic
  intermediate sequence grounded in the Politis/Romano/Wolf subsampling rate
  conditions, not a universally optimal finite-sample choice
- finite-sample results can be sensitive to `b`
- users should run sensitivity checks over a grid of `b`
- no method can certify validity from the data alone
- failed subsamples are diagnostic information, not something to hide

## Implementation Phases

### Phase 1: Core Generic Subsampling

- Add operation contract and private cache.
- Add ordinary observation-level without-replacement draw generator.
- Add public distribution, p-value, and CI methods.
- Support only `scaling = "sqrt_n"`.
- Support explicit integer `b` and deterministic `b = NULL`.
- Add simple mean-difference tests.

### Phase 2: Minimum-Volatility `b` Selection

- Add `select_optimal_b_subsampling()`.
- Adapt the minimum-volatility grid logic from
  `../PTE/PTE/R/select_optimal_m_prop.R`.
- Add `compute_subsampling_sensitivity()`.
- Add selector tests before using selection inside inference methods.

### Phase 3: Design-Aware Units

- Add unit resolver.
- Add blocking, cluster, matching, and KK matched/reservoir draw support.
- Add draw-contract tests for each design family.

### Phase 4: Diagnostics And Debug API

- Wire subsampling diagnostics into `compute_pval_debug()` and
  `compute_ci_debug()`.
- Preserve optimizer diagnostics for failed refits in debug mode when available.
- Add finite-fraction and failure-reason summaries.
- Add `b_selection` diagnostics when minimum-volatility selection is used.

### Phase 5: Path Audit Integration

- Add path-audit columns for subsampling p-values and CIs.
- Keep cells light green until class-specific empirical stability is measured.
- Add low-estimability summaries for subsampling separately from bootstrap.

## Acceptance Criteria

The feature is complete when:

- ordinary-design subsampling methods work for simple mean difference
- explicit `b`, deterministic `b = NULL`, and minimum-volatility `b`
  selection are implemented and tested
- design-aware draw tests pass for blocking, clustering, and matching
- p-value and CI methods return typed non-estimable output on invalid/unstable
  paths
- debug mode exposes finite fractions and failure reasons
- debug mode records the `b` selection grid and chosen `b` when applicable
- no existing `compute_bootstrap_*()` output changes
- no scalar fast path performs subsampling or optimizer diagnostics unless a
  subsampling method is explicitly called
