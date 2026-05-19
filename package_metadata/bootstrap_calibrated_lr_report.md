# Bootstrap-Calibrated Likelihood-Ratio Report

## Scope

This report covers the package's **likelihood-backed inference paths**, meaning
the concrete classes that participate in the `get_likelihood_test_spec()`
surface used by `InferenceAsympLik` for likelihood-ratio tests and confidence
interval inversion.

The question is:

> How hard would it be to implement a parametric-bootstrap-calibrated LR test
> and the corresponding inverted confidence interval across the
> likelihood-backed inference paths?

I use three labels:

- **Easy**: the path is a realistic first target for a parametric
  bootstrap-calibrated LR test.
- **Borderline**: the path is plausible, but needs family-specific simulation or
  more bespoke likelihood plumbing.
- **Difficult**: the path is too mixture-, quadrature-, frailty-, copula-, or
  custom-hybrid-heavy, or not naturally generative enough, for a practical
  package-wide rollout.

This is an implementation audit, not a statement about theoretical possibility.

## Short Answer

Bootstrap-calibrated LR is a promising idea for this repo, and the class
hierarchy now has the right structural home for it.

What the package already has:

- `InferenceParamBootstrap` as a dedicated intermediate class for families that
  may support parametric bootstrap LR calibration
- bootstrap infrastructure and reusable bootstrap workers
- LR test plumbing via `get_likelihood_test_spec()`
- constrained null refits with warm-start and memoization support

What it does **not** yet have:

- actual null-simulation machinery in `InferenceParamBootstrap` (the current
  stub just returns `FALSE` from `supports_lik_ratio_param_bootstrap()`)
- a generic hook to **simulate new responses under the fitted null likelihood**
- reusable worker support for **parametric null simulation**, not just
  resampling observed rows

So:

- the **structural foundation** in `InferenceParamBootstrap` is already in place
- implementing it across the **simple generative likelihood families** is
  realistic, since they all inherit via `InferenceAsympLikStdModCache`
- implementing it across **all** current LR paths is **hard** and intentionally
  out of scope for `InferenceParamBootstrap`

## Architecture

The inheritance chain relevant to parametric bootstrap is:

```
InferenceAsympLik
├── InferenceParamBootstrap            ← bootstrap-capable families live here
│   └── InferenceAsympLikStdModCache   ← GLM / KM standard-cache families
│       └── (Easy / Borderline concrete classes)
└── (direct InferenceAsympLik children) ← Difficult families; bypass InferenceParamBootstrap
```

The **Difficult** families (Cox, stratified Cox, combined-design Cox, frailty,
copula, clogit-plus-GLMM hybrid, custom combined-likelihood paths) remain direct
children of `InferenceAsympLik` and never pass through
`InferenceParamBootstrap`. This is by design: the intermediate class only exists
where parametric null simulation is a realistic long-run goal.

All **Easy** and most **Borderline** families already sit under
`InferenceParamBootstrap` via `InferenceAsympLikStdModCache`.

## What A Bootstrap-Calibrated LR Test Needs

For a null value `delta`, the observed LR statistic is:

```text
LR_obs(delta) = 2 { l(theta_hat) - l(theta_hat_delta) }
```

where:

- `theta_hat` is the unrestricted fit
- `theta_hat_delta` is the constrained fit under the null treatment effect

A bootstrap-calibrated LR test then:

1. fits the null model at `delta`
2. simulates bootstrap datasets from that fitted null model
3. recomputes the same LR statistic on each bootstrap dataset
4. estimates the p-value by the empirical tail probability

That is different from the package's current bootstrap workflows, which are
primarily:

- case/design resampling of observed rows
- bootstrap worker reuse for observed-data resampling

Those are useful building blocks, but they are not yet a parametric null
bootstrap.

## What Can Be Centralized In InferenceParamBootstrap

The shared outer logic can live in `InferenceParamBootstrap`.

`InferenceParamBootstrap` inherits `InferenceAsympLik` and already knows how to:

- obtain `spec = get_likelihood_test_spec()`
- evaluate the observed full and null neg-log-likelihoods
- refit the constrained null model at a chosen `delta`
- cache and warm-start repeated null fits

So a generic method is plausible:

```r
compute_lik_ratio_bootstrap_two_sided_pval = function(delta = 0, B = 199, ...)
```

The generic part should:

1. compute the observed LR statistic
2. fit the observed-data null model at `delta`
3. ask the subclass to simulate one bootstrap dataset under that null fit
4. refit unrestricted and null models on the bootstrap dataset
5. compute `LR*`
6. return the empirical tail probability

The missing piece is a family-specific hook to be defined on
`InferenceParamBootstrap` (with default `NULL` return):

```r
simulate_under_lik_null = function(spec, delta, null_fit, worker_state){
  NULL
}
```

That hook must create a new response vector, and for some families also a new
`dead` indicator or latent mixture structure, while preserving the design matrix
and treatment assignments.

## Why The CI Is Much Harder Than The P-Value

A single bootstrap-calibrated LR p-value at one `delta` is expensive but
manageable.

A confidence interval by inversion is much more expensive, because for each
candidate `delta` visited during inversion you may need:

1. a constrained fit on the observed data
2. `B` null-simulated bootstrap datasets
3. an unrestricted and constrained refit on each bootstrap dataset

So if the p-value is feasible, the CI is conceptually straightforward but
computationally much heavier.

This means the report's difficulty labels should be read as:

- "easy enough to support bootstrap-LR p-values first"
- with CI inversion being the same method but at significantly higher cost

## Main Architectural Gap

The current bootstrap worker abstraction is built around **resampling observed
units**:

- `bootstrap_sample_indices(...)`
- `load_bootstrap_sample_into_worker(...)`
- `bootstrap_subset_inference(...)`

Bootstrap-calibrated LR needs a second branch of worker logic inside
`InferenceParamBootstrap`:

- duplicate the inference object / design structure once
- keep `X`, `w`, and any matching/blocking structure fixed
- replace the response under a null-data simulator
- recompute unrestricted and constrained likelihood fits

So the biggest missing abstraction is not LR refitting (already in place via
`InferenceAsympLik`). It is **null-data simulation under each likelihood
family**, which is the gap that `InferenceParamBootstrap` needs to fill.

## Audit Table

### Incidence

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceIncidLogRegr` | `fast_logistic_regression_cpp` | **Easy** | Null simulation is straightforward Bernoulli generation from the constrained fit. |
| `InferenceIncidProbitRegr` | `fast_ordinal_probit_regression_cpp` in the 2-category case | **Easy** | Same as above with probit probabilities. |
| `InferenceIncidModifiedPoisson` | `fast_poisson_regression_cpp` | **Easy** | The LR path already uses a simple Poisson engine, so null simulation is straightforward. |
| `InferenceIncidKKModifiedPoisson` | `fast_poisson_regression_cpp` | **Easy** | Same Poisson null simulator as the companion path. |
| `InferenceIncidLogBinomial` | `fast_log_binomial_regression_cpp` | **Easy** | Bernoulli generation under the constrained fit is still straightforward. |
| `InferenceIncidBinomialIdentityRiskDiff` | `fast_identity_binomial_regression_cpp` | **Easy** | Same comment: if the constrained fit returns valid probabilities, null simulation is direct. |
| `InferenceIncidKKClogitOneLik` | stacked conditional-logistic + reservoir-logistic path via `fast_logistic_regression_with_var_cpp` | **Borderline** | The combined likelihood is still generative enough to simulate, but the data-generation interpretation is custom. |
| `InferenceIncidKKGLMM` | `fast_logistic_glmm_cpp` | **Borderline** | Parametric simulation is possible, but requires drawing random effects and respecting the GH-approximated likelihood structure. |
| `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm_cpp` | **Difficult** | Hybrid conditional-logit plus GLMM null simulation is too bespoke for a clean first rollout. Remains a direct `InferenceAsympLik` child. |

### Count

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceCountPoisson` | `fast_poisson_regression_cpp` | **Easy** | Direct Poisson null generation. |
| `InferenceCountRobustPoisson` | `fast_poisson_regression_cpp` | **Easy** | The LR path is still ordinary Poisson. |
| `InferenceCountQuasiPoisson` | `fast_poisson_regression_cpp` | **Easy** | Same comment: LR path is Poisson even if reported variance is quasi. |
| `InferenceCountNegBin` | `fast_neg_bin_cpp`, `fast_neg_bin_with_var_cpp` | **Borderline** | Null simulation is feasible, but dispersion handling adds family-specific work. |
| `InferenceCountZeroInflatedPoisson` | `fast_zero_augmented_poisson_cpp` | **Borderline** | Simulation is possible from the fitted mixture, but requires careful handling of the inflation component. |
| `InferenceCountHurdlePoisson` | `fast_zero_augmented_poisson_cpp` | **Borderline** | Same issue with added hurdle/truncation logic. |
| `InferenceCountZeroInflatedNegBin` | `fast_zinb_cpp` | **Borderline** | Generative, but mixture plus dispersion makes the simulator and refits more bespoke. |
| `InferenceCountHurdleNegBin` | `fast_hurdle_negbin_cpp` | **Borderline** | Same as above. |
| `InferenceCountKKGLMM` | `fast_poisson_glmm_cpp` | **Borderline** | Parametric null simulation is plausible but random effects make it materially more complex. |
| `InferenceCountKKHurdlePoissonOneLik` | `fast_hurdle_poisson_glmm_cpp` | **Difficult** | Hurdle plus GLMM structure is too custom for a broad first implementation. Remains a direct `InferenceAsympLik` child. |
| `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined_with_var_cpp` | **Difficult** | Custom combined likelihood with unclear generic null-simulation semantics. Remains a direct `InferenceAsympLik` child. |

### Continuous

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceContinKKOLSOneLik` | `fast_ols_with_var_cpp` | **Easy** | Gaussian null simulation is straightforward if the constrained fit exposes residual variance. |
| `InferenceContinKKGLMM` | `fast_gaussian_lmm_cpp` | **Borderline** | Parametric LMM simulation is feasible, but variance components and random effects add work. |

### Ordinal

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceOrdinalPropOddsRegr` | `fast_ordinal_regression_cpp` | **Borderline** | Implemented via ordinal category sampling from the constrained null fit. |
| `InferenceOrdinalOrderedProbitRegr` | `fast_ordinal_probit_regression_cpp` | **Borderline** | Implemented via ordinal category sampling from probit cumulative probabilities. |
| `InferenceOrdinalCauchitRegr` | `fast_ordinal_cauchit_regression_cpp` | **Borderline** | Implemented via ordinal category sampling from cauchit cumulative probabilities. |
| `InferenceOrdinalCloglogRegr` | `fast_ordinal_cloglog_regression_cpp` | **Borderline** | Implemented via ordinal category sampling from cloglog cumulative probabilities. |
| `InferenceOrdinalKKGLMM` | `fast_ordinal_glmm_cpp` | **Difficult** | Ordinal thresholds plus random effects move this out of the easy bootstrap-LR tier. Remains a direct `InferenceAsympLik` child. |

### Proportion

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferencePropBetaRegr` | `fast_beta_regression_cpp` | **Borderline** | Beta simulation is feasible, but the mean/precision parameterization needs family-specific handling. |
| `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta_cpp` | **Borderline** | Still generative, but three-part mixture simulation is more bespoke. |
| `InferencePropKKGLMM` | `fast_logistic_glmm_cpp` | **Borderline** | Same random-effects simulation issue as the other GLMM families. |

### Survival

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceSurvivalCoxPHRegr` | `fast_coxph_regression_cpp` | **Difficult** | The likelihood path is partial likelihood, not a full generative survival model. Null simulation requires extra baseline-hazard and censoring choices. Remains a direct `InferenceAsympLik` child. |
| `InferenceSurvivalStratCoxPHRegr` | `fast_stratified_coxph_regression_cpp` | **Difficult** | Same Cox issue with added strata-specific baseline structure. Remains a direct `InferenceAsympLik` child. |
| `InferenceSurvivalKKLWACoxOneLik` | `fast_coxph_regression_cpp` | **Difficult** | Combined-design Cox path still lacks a natural generic full-data simulator. Remains a direct `InferenceAsympLik` child. |
| `InferenceSurvivalKKStratCoxOneLik` | `fast_stratified_coxph_regression_cpp` | **Difficult** | Same issue with stratification. Remains a direct `InferenceAsympLik` child. |
| `InferenceSurvivalWeibullRegr` | `fast_weibull_regression_cpp` | **Borderline** | Implemented via Weibull null simulation with observed censoring carried forward. |
| `InferenceSurvivalDepCensTransformRegr` | `fast_dep_cens_transform_optim_cpp` | **Difficult** | Highly bespoke joint event/censoring structure. Remains a direct `InferenceAsympLik` child. |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty_cpp` | **Difficult** | Frailty simulation plus censoring plus repeated null refits is too heavy for a first generic rollout. Remains a direct `InferenceAsympLik` child. |
| `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_clayton_weibull_aft_optim_cpp` | **Difficult** | Copula-based survival simulation is possible in principle but too bespoke package-wide. Remains a direct `InferenceAsympLik` child. |

## Easy Tier

These are the best first targets:

1. `InferenceIncidLogRegr`
2. `InferenceIncidProbitRegr`
3. `InferenceCountPoisson`
4. `InferenceIncidModifiedPoisson`
5. `InferenceIncidKKModifiedPoisson`
6. `InferenceIncidLogBinomial`
7. `InferenceIncidBinomialIdentityRiskDiff`
8. `InferenceCountRobustPoisson`
9. `InferenceCountQuasiPoisson`
10. `InferenceContinKKOLSOneLik`

Why these are easiest:

- the null model is clearly generative
- simulating responses under the constrained fit is straightforward
- the unconstrained and constrained refits already exist
- all ten already sit under `InferenceParamBootstrap` (via
  `InferenceAsympLikStdModCache`), so the hook can be added without any
  class-hierarchy surgery

## Borderline Tier

These are plausible, but need family-specific simulators and more careful
worker-state support:

- `InferenceIncidKKClogitOneLik`
- `InferenceIncidKKGLMM`
- `InferenceCountNegBin`
- zero-inflated / hurdle count families
- `InferenceCountKKGLMM`
- `InferenceContinKKGLMM`
- ordinal fixed-effects models
- beta and zero/one-inflated beta
- `InferencePropKKGLMM`
- `InferenceSurvivalWeibullRegr`

The recurring complications are:

- dispersion / precision parameters
- mixture or hurdle generation
- random effects
- ordinal category sampling
- event-time and censoring generation

These are good second-tier candidates if the base framework proves stable.
All of them already inherit via `InferenceParamBootstrap`, so no hierarchy
changes are needed, only the null-simulation hook implementation.

## Difficult Tier

These are the weakest candidates for a common bootstrap-LR rollout and
**remain direct children of `InferenceAsympLik`**, bypassing
`InferenceParamBootstrap` entirely:

- Cox / stratified Cox / combined-design Cox
- the dependent-censoring transform path
- frailty survival models
- copula survival models
- hurdle-GLMM and other custom hybrid combined-likelihood paths

The recurring issue is that the existing likelihood test does not come with a
simple generic full-data simulator under the null.

In particular, the Cox paths are difficult because partial likelihood is not by
itself a complete generative model for survival times and censoring.

## Implementation Status

### Phase 1: InferenceParamBootstrap API — **DONE**

- `compute_lik_ratio_bootstrap_two_sided_pval(delta, B, ...)` implemented with
  parallel-chunk strategy (same multicore approach as `compute_randomization_two_sided_pval`)
- Hooks `supports_lik_ratio_param_bootstrap()` and `simulate_under_lik_null(spec, delta, null_fit)`
  added to `InferenceParamBootstrap`

### Phase 2: Easy Generative Families — **DONE**

All ten Easy families have `simulate_under_lik_null` and
`supports_lik_ratio_param_bootstrap = TRUE`:
`InferenceIncidLogRegr`, `InferenceIncidProbitRegr`, `InferenceCountPoisson`,
`InferenceIncidModifiedPoisson`, `InferenceIncidKKModifiedPoisson`,
`InferenceIncidLogBinomial`, `InferenceIncidBinomialIdentityRiskDiff`,
`InferenceCountRobustPoisson`, `InferenceCountQuasiPoisson`,
`InferenceContinKKOLSOneLik`.

### Phase 3: Selected Borderline Families — **DONE** (with noted exceptions)

| Family | Status | Notes |
|---|---|---|
| `InferenceCountNegBin` | ✓ Done | dispersion `theta` extracted from `null_fit$theta_hat` |
| `InferenceCountZeroInflatedPoisson` | ✓ Done | mixture simulated via `rbinom`/`rpois`; fixed pre-existing `neg_loglik` extraction bug |
| `InferenceCountHurdlePoisson` | ✓ Done | zero-truncated Poisson generation via `qpois` |
| `InferenceCountZeroInflatedNegBin` | ✓ Done | `log_theta` extracted from tail of params vector |
| `InferenceCountHurdleNegBin` | ✗ Disabled | C-level heap corruption in `fast_truncated_negbin_count_cpp` during loop; `supports_lik_ratio_param_bootstrap = FALSE` |
| `InferenceCountKKGLMM` | ✓ Done | random effects drawn before Poisson generation; inherits from `InferenceParamBootstrap` |
| `InferenceContinKKGLMM` | ✓ Done | all params in `$b = c(betas, log_sigma, log_tau)`; inherits from `InferenceParamBootstrap` |
| `InferencePropZeroOneInflatedBetaRegr` | ✓ Done | three-part mixture (zero / beta / one) simulation |
| `InferenceIncidKKClogitOneLik` | ✓ Done | generative Bernoulli simulation from combined-likelihood design matrix; fixed pre-existing `attempt$X_fit` → `attempt$X` bug |
| `InferenceIncidKKGLMM` | ✗ Skipped | alias for `InferenceIncidKKClogitPlusGLMMOneLik` (Difficult) |
| `InferencePropKKGLMM` | ✗ Skipped | `InferenceAbstractKKLogisticGLMMOneLik` — combined clogit+GLMM (Difficult) |
| Ordinal fixed-effects models | ✓ Done | `InferenceOrdinalPropOddsRegr`, `InferenceOrdinalOrderedProbitRegr`, `InferenceOrdinalCauchitRegr`, and `InferenceOrdinalCloglogRegr` now simulate ordinal categories under the constrained null fit |
| `InferenceSurvivalWeibullRegr` | ✓ Done | Weibull event times simulated under the constrained null fit with observed censoring thresholds carried forward |

## Concrete Checklist

### Checklist A: Finish The Parametric-Bootstrap LR Architecture

- [x] Add the generic outer-loop API in `InferenceParamBootstrap`:
  `compute_lik_ratio_bootstrap_two_sided_pval(...)`
- [x] Add bootstrap-LR CI inversion in `InferenceParamBootstrap`:
  `compute_lik_ratio_bootstrap_confidence_interval(...)`
- [x] Add the family hook `simulate_under_lik_null(spec, delta, null_fit)`
- [x] Add the family capability gate `supports_lik_ratio_param_bootstrap()`
- [x] Add a validator for the return object from `simulate_under_lik_null(...)`
  so all families return the same minimal contract:
  `full_fit`, `fit_null(delta, start = NULL)`, and `neg_loglik(fit)`
- [x] Add reusable simulation helpers inside `InferenceParamBootstrap` or a
  nearby helper module for common null generators:
  Bernoulli, Poisson, Gaussian, ordinal-category sampling, Weibull
- [x] Add explicit replicate diagnostics for bootstrap-LR failures:
  simulated-data failure, full-refit failure, null-refit failure, non-finite LR
- [x] Add a retry / minimum-usable-replicates policy when too many bootstrap
  replicates return `NA`
- [ ] Audit every `InferenceParamBootstrap` descendant and mark each one as one
  of:
  supported now, intentionally unsupported, or disabled pending backend fixes
- [ ] Keep `InferenceCountHurdleNegBin` disabled until the
  `fast_truncated_negbin_count_cpp` stability issue is fixed

### Checklist B: Add Reusable Worker Support For Parametric Null Simulation

- [x] Add a parametric-bootstrap worker-state API parallel to the existing
  nonparametric worker API, likely in `InferenceParamBootstrap`:
  `create_param_bootstrap_worker_state(...)`
- [ ] Add a worker loader for simulated null data:
  `load_param_bootstrap_draw_into_worker(worker_state, sim_data)`
- [x] Add a worker computation method:
  `compute_param_bootstrap_worker_lrt(worker_state, delta)`
- [ ] Make worker state hold fixed design-side objects once:
  `X`, `w`, block / match structure, static metadata, base warm starts
- [ ] On each bootstrap replicate, replace only response-side data:
  `y`, `dead`, and any family-specific latent / mixture components
- [ ] Refit unrestricted and null models on the worker object using the same
  family-specific likelihood machinery already used in the generic path
- [x] Route `compute_lik_ratio_bootstrap_two_sided_pval(...)` through the worker
  path when available, with fallback to the current per-replicate path
- [x] Add parity tests:
  reusable-worker path vs current generic path for the same seed and `B`
- [ ] Add serial / parallel determinism tests where exact equality is expected
  under fixed seeds
- [ ] Add smoke tests for censoring, mixture, ordinal, and GLMM-style families

### Checklist C: Make The Existing User-Facing API Easy To Run Reliably

- [ ] Keep the existing public methods as the sole user-facing API:
  `compute_lik_ratio_bootstrap_two_sided_pval(...)` and
  `compute_lik_ratio_bootstrap_confidence_interval(...)`
- [ ] Standardize and document the expected user-facing arguments and defaults:
  `delta = 0`, `B = 199`, `show_progress = FALSE`
- [ ] Add roxygen docs describing:
  what each method does, when it is available, and expected runtime costs
- [ ] Document the availability rule clearly:
  only classes with `supports_lik_ratio_param_bootstrap() == TRUE` should
  support successful runs
- [ ] Add examples to a few stable first-wave families:
  `InferenceIncidLogRegr`, `InferenceCountPoisson`,
  `InferenceSurvivalWeibullRegr`, `InferenceOrdinalPropOddsRegr`
- [ ] Add one package-level vignette / README section showing the end-user flow:
  fit design -> create inference object -> call bootstrap LR p-value
- [ ] Decide whether unsupported families should:
  error, return `NA`, or advertise non-support via a helper
- [ ] Optionally add a helper such as
  `supports_bootstrap_lrt()` for easy user-side capability checks

### Checklist D: What A User Needs In Order To Run The Existing API Reliably

- [ ] Pick a first-wave supported family
- [ ] Fit the corresponding inference object on a completed design
- [ ] Confirm the class is one with `supports_lik_ratio_param_bootstrap() == TRUE`
- [ ] Start with a small `B` smoke run, e.g. `B = 19` or `B = 49`
- [ ] Then run the real analysis at a larger `B`, e.g. `B = 199` or `B = 499`
- [ ] For expensive families, prefer p-values first and only then CI inversion
- [ ] Expect CI inversion to be much slower than the p-value because each
  candidate `delta` triggers a full nested bootstrap-LR calculation
- [ ] For families with occasional optimizer failures, monitor the share of
  usable bootstrap replicates and rerun with a larger `B` if needed

## Bottom Line

Bootstrap-calibrated LR testing is a realistic extension for this repo.
The `InferenceParamBootstrap` class is now the designated home for this
functionality, and the class hierarchy already encodes the easy/difficult split:

- **Easy / Borderline** families sit under `InferenceParamBootstrap` and are
  ready to receive the null-simulation hook
- **Difficult** families remain direct children of `InferenceAsympLik` and are
  explicitly excluded

The actual gap is no longer structural. The generic parametric-null bootstrap
loop in `InferenceParamBootstrap` and the main first-wave
`simulate_under_lik_null(...)` implementations are now in place. What remains
is selective hardening, test coverage, and deliberate scoping decisions for the
most bespoke hybrid families.
