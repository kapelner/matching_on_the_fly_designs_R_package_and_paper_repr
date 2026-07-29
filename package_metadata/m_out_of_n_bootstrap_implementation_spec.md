# m-out-of-n Bootstrap Implementation Spec

Generated: 2026-07-27

## Scope

This spec defines an implementation plan for m-out-of-n bootstrap inference in
EDI.

The m-out-of-n bootstrap draws `m < n` exchangeable units with replacement,
recomputes the estimator on each size-`m` bootstrap sample, and uses a
rate-corrected distribution to approximate the full-sample estimator's sampling
law.

This is not the existing ordinary nonparametric bootstrap. The current
bootstrap infrastructure draws `n` units with replacement. The m-out-of-n
bootstrap keeps the with-replacement bootstrap mechanism but changes the
resample size and calibration.

Related implementation documents:

- [prw_subsampling_implementation_spec.md](prw_subsampling_implementation_spec.md)
- [public_diagnostics_api_spec.md](public_diagnostics_api_spec.md)
- [optimizer_diagnostics_report.md](optimizer_diagnostics_report.md)
- [../package_tests/path_audit_hardening_report.md](../package_tests/path_audit_hardening_report.md)

## Motivation

m-out-of-n bootstrap gives EDI a principled alternative when the ordinary
n-out-of-n bootstrap is unstable or theoretically questionable for nonregular,
boundary-sensitive, or nonsmooth estimators.

It should be exposed as an explicit inference family. It should not silently
replace ordinary bootstrap when ordinary bootstrap has failed.

## Non-Goals

- Do not change existing `compute_bootstrap_*()` behavior.
- Do not silently drop failed ordinary bootstrap draws and call the result
  m-out-of-n inference.
- Do not implement PRW subsampling in this spec. PRW subsampling draws without
  replacement and has its own API and operation contract.
- Do not claim universal validity. m-out-of-n bootstrap still requires a
  defensible choice of `m` and method-specific regularity assumptions.

## Public API

Add explicit public methods:

```r
approximate_m_out_of_n_bootstrap_distribution_beta_hat_T(
  B = 501,
  m = NULL,
  show_progress = TRUE,
  debug = FALSE,
  bootstrap_type = NULL,
  scaling = "sqrt_n",
  center = "full_estimate"
)

compute_m_out_of_n_bootstrap_two_sided_pval(
  delta = 0,
  B = 501,
  m = NULL,
  type = "centered",
  show_progress = TRUE,
  min_number_usable_samples = 5L,
  bootstrap_type = NULL,
  scaling = "sqrt_n"
)

compute_m_out_of_n_bootstrap_confidence_interval(
  alpha = 0.05,
  B = 501,
  m = NULL,
  type = "basic",
  show_progress = TRUE,
  min_number_usable_samples = 5L,
  bootstrap_type = NULL,
  scaling = "sqrt_n"
)
```

Do not add `m` to `compute_bootstrap_two_sided_pval()` or
`compute_bootstrap_confidence_interval()` in the first pass. The tuning and
calibration are different enough that users should opt into the method by
name.

## Theory Contract

For a full-sample estimate `theta_hat_n`, an m-out-of-n bootstrap estimate
`theta_hat_m_star`, and a scaling sequence `tau_n`, estimate the distribution
of:

```text
tau_n * (theta_hat_n - theta)
```

with the empirical distribution of:

```text
tau_m * (theta_hat_m_star - theta_hat_n)
```

where `m -> infinity` and, in nonregular settings, typically `m / n -> 0`.

For root-n regular estimators, the default scaling is:

```text
tau_n = sqrt(n)
tau_m = sqrt(m)
```

The confidence interval should use quantiles of:

```text
T_m_star = tau_m * (theta_hat_m_star - theta_hat_n)
```

and return:

```text
CI = c(theta_hat_n - q_hi / tau_n,
       theta_hat_n - q_lo / tau_n)
```

The centered two-sided p-value for `H0: theta = delta` should compare:

```text
T_n = tau_n * (theta_hat_n - delta)
```

to `T_m_star` using the same empirical two-sided tail convention as PRW
subsampling.

## Shared Minimum-Volatility Selector

The minimum-volatility approach from PTE
(`../PTE/PTE/R/select_optimal_m_prop.R`) should be reused conceptually for both
m-out-of-n bootstrap and PRW subsampling.

The shared idea is:

1. evaluate a grid of candidate size exponents
2. compute an inference answer at each grid point
3. measure local volatility over a complete centered neighborhood
4. select the eligible interior grid point whose answer is most stable
5. rerun final inference at the selected size with production `B`

For m-out-of-n bootstrap, the tuning size is:

```r
m = floor(n_units ^ gamma)
```

For PRW subsampling, the same selector engine uses:

```r
b = floor(n_units ^ gamma)
```

The selector engine should be generic:

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

The selector logic is shared; the resampling draw and calibration are not.
m-out-of-n bootstrap draws with replacement, while PRW subsampling draws
without replacement.

## Choosing `m`

The first implementation should accept:

```r
m = <integer>
m = NULL
m = list(method = "minimum_volatility", ...)
```

The deterministic fallback rule for `m = NULL` is:

```r
m = floor(n_units ^ 0.7)
```

with guards:

```r
m >= max(5L, p_effective + 2L)
m <= floor(n_units / 2)
```

For production tuning, add:

```r
select_optimal_m_out_of_n_bootstrap(
  B = 251,
  alpha = 0.05,
  m_pow_of_n_grid = seq(0.5, 0.9, by = 0.05),
  m_grid = NULL,
  objective = "ci_width",
  target = "ci",
  volatility_window = 3L,
  bootstrap_type = NULL,
  scaling = "sqrt_n",
  show_progress = TRUE,
  min_finite_fraction = 0.8
)
```

If `m_grid` is supplied, use it directly. Otherwise compute:

```r
m_grid = unique(floor(n_units ^ m_pow_of_n_grid))
```

Grid points that fail validity checks should remain in `grid_table` with a
typed failure reason, but they are not eligible for selection.

The first selector implementation should support:

```r
objective = "ci_width"
target = "ci"
```

Additional objectives can be added later:

```r
objective = "pval_stability"
objective = "estimate_stability"
objective = "finite_fraction_penalized_ci_width"
```

For the default objective, define volatility as:

```r
sd(ci_width[(k - h):(k + h)], na.rm = FALSE)
```

where `h = (volatility_window - 1) / 2`. Boundary grid points are not eligible,
matching the PTE implementation.

Return:

```r
list(
  m_optimal = <integer>,
  m_pow_of_n_optimal = <numeric or NA>,
  grid_table = <data.frame>,
  objective = <character>,
  target = <character>,
  B = <integer>,
  alpha = <numeric>,
  volatility_window = <integer>,
  bootstrap_type = <character or NULL>,
  scaling = <character>
)
```

with class:

```r
"EDIMOutOfNBootstrapMSelection"
```

For `m = NULL`, use the deterministic rule so ordinary calls remain cheap.
For:

```r
m = list(method = "minimum_volatility")
```

run the selector first, then rerun final inference at `m_optimal`. Do not reuse
the low-`B` calibration grid as final reported inference unless the user
explicitly requested the same `B`, `m`, and arguments.

## Design-Aware Resampling Units

Use the same exchangeable-unit resolver planned for PRW subsampling. The
difference is replacement:

- m-out-of-n bootstrap samples `m` units with replacement
- PRW subsampling samples `b` units without replacement

For clusters, pairs, matched sets, and blocks, the draw may contain multiple
rows per selected unit and duplicate units are allowed because this is still a
bootstrap.

For blocking designs, preserve the existing `bootstrap_type` vocabulary where
possible:

```r
bootstrap_type = NULL
bootstrap_type = "within_blocks"
bootstrap_type = "resample_blocks"
```

but make the size contract explicit in terms of exchangeable units.

## Internal Architecture

Add a new operation contract:

```r
m_out_of_n_boot = list(
  operation = "m_out_of_n_boot",
  draw_type = "unit_sample_with_replacement_size_m",
  loader = "load_m_out_of_n_bootstrap_draw_into_worker",
  estimator = "compute_bootstrap_worker_estimate",
  cache_name = "m_out_of_n_boot_distr_cache",
  cache_key_method = "m_out_of_n_bootstrap_cache_key"
)
```

Do not reuse `"non_param_boot"` because cache keys and calibration differ.

## Draw Representation

Represent each draw as:

```r
list(
  i_b = <integer row indices>,
  m_vec_b = <integer or NULL>,
  unit_ids = <integer unit ids>,
  unit_type = <character>,
  m = <integer>,
  n_units = <integer>,
  stratum_ids = <integer or NULL>
)
```

The row index vector may contain duplicate exchangeable units. This is the
main draw-level distinction from PRW subsampling.

## Failure Handling

Required typed reasons:

```text
m_out_of_n_original_estimate_unavailable
m_out_of_n_invalid_resample_size
m_out_of_n_too_few_finite_estimates
m_out_of_n_unit_structure_unavailable
m_out_of_n_extreme_confidence_interval
m_out_of_n_m_selection_failed
m_out_of_n_m_selection_too_few_finite_estimates
```

The method may remove non-finite replicate estimates only if the finite count
and finite fraction pass explicit thresholds. The debug API should expose
failure reasons and finite fractions.

## Diagnostics Integration

`compute_pval_debug()` and `compute_ci_debug()` should expose:

```r
diagnostics$m_out_of_n_bootstrap
```

Required fields:

```r
list(
  n_units = <integer>,
  m = <integer>,
  B = <integer>,
  scaling = <character>,
  center = "full_estimate",
  unit_type = <character>,
  n_attempted = <integer>,
  n_finite = <integer>,
  finite_fraction = <numeric>,
  failure_reasons = <table/list>,
  q_lo = <numeric or NA>,
  q_hi = <numeric or NA>,
  m_selection = <list or NULL>
)
```

## Tests

Add focused tests before path-audit integration:

- ordinary draws contain exactly `m` sampled units with replacement
- duplicate units are possible and represented correctly
- cluster/pair/matched-set draws preserve intact unit rows
- explicit `m` is respected
- deterministic `m = NULL` is bounded correctly
- `m_pow_of_n_grid` maps to unique integer `m` values
- boundary grid points are ineligible for minimum-volatility selection
- selector chooses the least volatile eligible interior point
- tie handling chooses the smallest `m`
- failed grid points remain in `grid_table` with typed reasons
- final inference with `m = list(method = "minimum_volatility")` reruns at
  `m_optimal`

Use simple mean difference as the first inference target, then add one benign
model-fit path and one sparse path with failed resamples.

## Documentation

Document m-out-of-n bootstrap separately from ordinary bootstrap and PRW
subsampling.

Required caveats:

- validity depends on `m -> infinity` and, in nonregular cases, `m / n -> 0`
  following the m-out-of-n bootstrap framework of Bickel, Gotze, and van Zwet
  (1997) and the size-selection discussion in Bickel and Sakov (2008)
- the default `m = NULL` uses `floor(n_units ^ 0.7)`, a deterministic
  intermediate sequence satisfying the asymptotic rate conditions, not a
  universally optimal finite-sample choice
- finite-sample results can be sensitive to `m`
- minimum-volatility selection is a calibration heuristic, not a proof of
  validity
- users should inspect the selection grid and finite fractions
- failed replicates are diagnostic information

## Implementation Phases

### Phase 1: Core m-out-of-n Bootstrap

- Add operation contract and private cache.
- Add ordinary unit-level with-replacement size-`m` draw generator.
- Add public distribution, p-value, and CI methods.
- Support explicit integer `m` and deterministic `m = NULL`.
- Support only `scaling = "sqrt_n"`.

### Phase 2: Shared Minimum-Volatility Selector

- Add generic `select_optimal_resample_size()`.
- Add `select_optimal_m_out_of_n_bootstrap()`.
- Reuse the same selector engine later from
  `select_optimal_b_subsampling()`.
- Add selector tests.

### Phase 3: Design-Aware Units

- Add cluster, blocking, matching, and KK matched/reservoir draw support.
- Keep replacement semantics explicit in tests.

### Phase 4: Diagnostics And Debug API

- Wire m-out-of-n diagnostics into the public debug API.
- Preserve optimizer diagnostics for failed refits in debug mode when
  available.

### Phase 5: Path Audit Integration

- Add separate path-audit columns for m-out-of-n p-values and CIs.
- Do not fold m-out-of-n into ordinary bootstrap columns.

## Acceptance Criteria

The feature is complete when:

- ordinary-design m-out-of-n methods work for simple mean difference
- explicit `m`, deterministic `m = NULL`, and minimum-volatility `m`
  selection are implemented and tested
- design-aware draw tests pass for blocking, clustering, and matching
- p-value and CI methods return typed non-estimable output on invalid/unstable
  paths
- debug mode records the `m` selection grid and chosen `m` when applicable
- no existing `compute_bootstrap_*()` output changes
