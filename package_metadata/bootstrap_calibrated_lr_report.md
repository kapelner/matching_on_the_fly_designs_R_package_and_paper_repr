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

## Current `InferenceParamBootstrap` Descendant Status Audit

Status labels in this section mean:

- **Supported now**: the class currently opts into
  `supports_lik_ratio_param_bootstrap()`, either unconditionally or when the
  Rcpp backend is enabled.
- **Intentionally unsupported**: the class is intentionally kept off the
  `InferenceParamBootstrap` branch, so bootstrap-LR is outside its current
  supported surface and the public bootstrap-LR methods are not exposed.
- **Disabled pending backend fixes**: the class has a bootstrap simulator path
  but is explicitly switched off because the backend is not currently stable
  enough to expose.

Abstract / infrastructure-only bootstrap descendants omitted from the concrete
audit: `InferenceAsympLikStdModCache`, `InferenceCountLikelihood`,
`InferenceCountCompositeLikelihood`, `InferenceCountZeroAugmentedPoissonAbstract`,
`InferenceKKPassThroughCompound`, `InferenceAbstractKKLWACoxOneLik`, and
`InferenceAbstractKKWeibullFrailtyOneLik`.

### Likelihood-Backed Concrete Descendants

| Class | Status | Notes |
|---|---|---|
| `InferenceIncidLogRegr` | Supported now | Bernoulli null simulation implemented. |
| `InferenceIncidProbitRegr` | Supported now | Bernoulli null simulation implemented via probit link. |
| `InferenceIncidModifiedPoisson` | Supported now | Poisson null simulation implemented. |
| `InferenceIncidLogBinomial` | Supported now | Bernoulli null simulation implemented. |
| `InferenceIncidBinomialIdentityRiskDiff` | Supported now | Bernoulli null simulation implemented with identity-link clipping checks. |
| `InferenceIncidKKCondLogitOneLik` | Supported now | Combined-likelihood Bernoulli simulation implemented. |
| `InferenceCountPoisson` | Supported now | Poisson null simulation implemented. |
| `InferenceCountRobustPoisson` | Supported now | Inherits the supported `InferenceCountCompositeLikelihood` bootstrap path. |
| `InferenceCountQuasiPoisson` | Supported now | Inherits the supported `InferenceCountCompositeLikelihood` bootstrap path. |
| `InferenceCountNegBin` | Supported now | Negative-binomial null simulation implemented. |
| `InferenceCountZeroInflatedPoisson` | Supported now | Zero-inflated Poisson simulation implemented. |
| `InferenceCountHurdlePoisson` | Supported now | Hurdle-Poisson simulation implemented. |
| `InferenceCountZeroInflatedNegBin` | Supported now | Zero-inflated negative-binomial simulation implemented. |
| `InferenceCountKKGLMM` | Supported now | Rcpp-backed only; `supports_lik_ratio_param_bootstrap()` is gated by `use_rcpp`. |
| `InferenceContinKKOLSOneLik` | Supported now | Gaussian null simulation implemented. |
| `InferenceContinKKGLMM` | Supported now | Rcpp-backed only; `supports_lik_ratio_param_bootstrap()` is gated by `use_rcpp`. |
| `InferenceOrdinalPropOddsRegr` | Supported now | Ordinal-category simulation implemented. |
| `InferenceOrdinalOrderedProbitRegr` | Supported now | Ordinal-category simulation implemented. |
| `InferenceOrdinalCauchitRegr` | Supported now | Ordinal-category simulation implemented. |
| `InferenceOrdinalCloglogRegr` | Supported now | Ordinal-category simulation implemented. |
| `InferencePropBetaRegr` | Supported now | Beta null simulation implemented. |
| `InferencePropZeroOneInflatedBetaRegr` | Supported now | Zero/one/beta mixture simulation implemented. |
| `InferenceSurvivalCoxPHRegr` | Supported now | Rcpp-backed only; baseline-hazard-based simulator implemented behind `use_rcpp`. |
| `InferenceSurvivalStratCoxPHRegr` | Supported now | Rcpp-backed only; stratified baseline-hazard simulator implemented behind `use_rcpp`. |
| `InferenceSurvivalWeibullRegr` | Supported now | Weibull null simulation implemented. |
| `InferenceSurvivalKKLWACoxPHOneLik` | Supported now | Combined-design Cox null simulation implemented. |
| `InferenceSurvivalKKStratCoxPHOneLik` | Supported now | Combined-design stratified Cox null simulation implemented. |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | Supported now | Rcpp-backed only; frailty simulation implemented behind `use_rcpp`. |
| `InferenceSurvivalKKClaytonCopulaOneLik` | Supported now | Clayton-copula survival simulation implemented. |
| `InferenceCountHurdleNegBin` | Disabled pending backend fixes | Bootstrap simulator code exists, but the class is kept off `InferenceParamBootstrap` until the repeated null-refit path is reliable. |
| `InferenceIncidRiskDiff` | Intentionally unsupported | Kept off `InferenceParamBootstrap`; no likelihood-null simulator hook. |
| `InferencePropFractionalLogit` | Intentionally unsupported | Kept off `InferenceParamBootstrap`; not yet given a generative fractional-response bootstrap path. |
| `InferenceOrdinalAdjCatLogitRegr` | Intentionally unsupported | Kept off `InferenceParamBootstrap`; no adjacent-category null simulator yet. |
| `InferenceOrdinalStereotypeLogitRegr` | Intentionally unsupported | Kept off `InferenceParamBootstrap`. |
| `InferenceOrdinalContRatioRegr` | Intentionally unsupported | Kept off `InferenceParamBootstrap`. |
| `InferenceSurvivalDepCensTransformRegr` | Intentionally unsupported | Kept off `InferenceParamBootstrap`; the joint event/censoring simulator remains out of scope. |

### Pass-Through Concrete Descendants

These classes use KK pass-through layers that are intentionally kept off
`InferenceParamBootstrap`, so they do not expose the likelihood-ratio
parametric bootstrap methods.

| Class | Status | Notes |
|---|---|---|
| `InferenceAllKKMeanDiffIVWC` | Intentionally unsupported | Pass-through IVWC mean-difference path; no LR bootstrap hook. |
| `InferenceContinKKRobustRegrOneLik` | Intentionally unsupported | Robust one-likelihood path is kept off `InferenceParamBootstrap`. |
| `InferenceContinKKRobustRegrIVWC` | Intentionally unsupported | Pass-through IVWC robust regression path; no LR bootstrap hook. |
| `InferenceContinKKOLSIVWC` | Intentionally unsupported | IVWC OLS path does not expose a parametric LR bootstrap path. |
| `InferenceBaiAdjustedTKK14` | Intentionally unsupported | Bai-adjusted test path is not on the current LR-bootstrap surface. |
| `InferenceBaiAdjustedTKK21` | Intentionally unsupported | Bai-adjusted test path is not on the current LR-bootstrap surface. |
| `InferenceContinKKQuantileRegrOneLik` | Intentionally unsupported | Quantile-regression path does not currently define a likelihood-null simulator. |
| `InferenceContinKKQuantileRegrIVWC` | Intentionally unsupported | IVWC quantile path does not currently define a likelihood-null simulator. |
| `InferencePropKKQuantileRegrOneLik` | Intentionally unsupported | Proportion quantile path does not currently define a likelihood-null simulator. |
| `InferenceIncidKKNewcombeRiskDiff` | Intentionally unsupported | Newcombe-style risk-difference path is not on the LR-bootstrap surface. |

## Additional Families That Could Be Brought Onto The Parametric-Bootstrap Branch

### Best Candidates

- `InferenceOrdinalAdjCatLogitRegr`
  at `EDI/R/inference_ordinal_adj_cat_logit.R`
  already has a dedicated C++ fit path, warm starts, and reusable-worker
  support via `generate_mod()` and `get_bootstrap_worker_spec()`. What it
  lacks is the same likelihood-test spec plus constrained-null simulation hook
  already present in the supported ordinal families.
- `InferenceOrdinalStereotypeLogitRegr`
  at `EDI/R/inference_ordinal_stereotype_logit.R`
  has a real parametric backend, warm starts, and reusable-worker plumbing,
  but no LR-bootstrap contract yet.
- `InferenceOrdinalContRatioRegr`
  at `EDI/R/inference_ordinal_stereotype_logit.R`
  is in the same position as the stereotype model: proper parametric backend,
  no null-simulation/bootstrap-LR hook yet.
- `InferenceCountHurdleNegBin`
  at `EDI/R/inference_count_hurdle.R`
  is a deferred reactivation candidate rather than a fresh addition. It
  already has `get_likelihood_test_spec()` and `simulate_under_lik_null()`,
  but should stay off the branch until the repeated null-refit failure is
  fixed.

### Plausible But Harder

- `InferenceSurvivalDepCensTransformRegr`
  at `EDI/R/inference_survival_dep_cens_transform.R`
  is closer than most unsupported families because it already implements
  `supports_likelihood_tests()` and `get_likelihood_test_spec()`. The missing
  piece is `simulate_under_lik_null()`. Its backend model is explicit enough to
  simulate from: a joint lognormal event/censoring model with correlated
  normals and parameters
  `[beta_event, beta_cens, log_sigma_event, log_sigma_cens, atanh_rho]`.
  That makes it feasible, but materially trickier than the ordinal additions.

### Probably Not Good Candidates For This Branch

- `InferencePropFractionalLogit`
  at `EDI/R/inference_proportion_fractional_logit.R`
  explicitly sets `supports_likelihood_tests = FALSE`. This is
  quasi-likelihood territory, so there is no obvious native LR target without
  introducing a different generative family.
- `InferenceIncidRiskDiff`
  at `EDI/R/inference_incidence_risk_diff.R`
  is OLS / linear-probability-model based rather than likelihood based, so LR
  parametric bootstrap is not a natural fit.
- The KK pass-through / IVWC / Bai / quantile / Newcombe families are not
  single parametric likelihoods, so putting them on the LR parametric-bootstrap
  branch would force the wrong abstraction.

### Suggested Next-Wave Order

1. `InferenceOrdinalAdjCatLogitRegr`
2. `InferenceOrdinalStereotypeLogitRegr`
3. `InferenceOrdinalContRatioRegr`
4. Re-enable `InferenceCountHurdleNegBin` after the null-refit bug is fixed
5. `InferenceSurvivalDepCensTransformRegr`
