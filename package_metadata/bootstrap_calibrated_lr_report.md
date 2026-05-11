# Bootstrap-Calibrated Likelihood-Ratio Report

## Scope

This report covers the package’s **likelihood-backed inference paths**, meaning
the concrete classes that participate in the `get_likelihood_test_spec()`
surface used by `InferenceAsympLik` for likelihood-ratio tests and confidence
interval inversion.

The question is:

> How hard would it be to implement a
> `compute_lik_ratio_bootstrap_two_sided_pval(...)` method, and the
> corresponding inverted confidence interval, across the likelihood-backed
> inference paths?

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

Bootstrap-calibrated LR is a promising idea for this repo, but the current
bootstrap support only gets you part of the way there.

What the package already has:

- bootstrap infrastructure
- reusable bootstrap workers
- LR test plumbing via `get_likelihood_test_spec()`
- constrained null refits
- warm-start and memoization support

What it does **not** yet have:

- a generic hook to **simulate new responses under the fitted null likelihood**
- reusable worker support for **parametric null simulation**, not just
  resampling observed rows
- a clear generative strategy for partial-likelihood and custom hybrid paths

So:

- implementing the **base framework** in `InferenceAsympLik` is **moderate**
- implementing it across the **simple generative likelihood families** is
  realistic
- implementing it across **all** current LR paths is **hard**

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

That is different from the package’s current bootstrap workflows, which are
primarily:

- case/design resampling of observed rows
- bootstrap worker reuse for observed-data resampling

Those are useful building blocks, but they are not yet a parametric null
bootstrap.

## What Can Be Centralized

The shared outer logic can live in `InferenceAsympLik`.

The base class already knows how to:

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

The missing piece is a family-specific hook like:

```r
simulate_under_lik_null = function(spec, delta, null_fit, worker_state){
  NULL
}
```

or equivalently:

```r
generate_bootstrap_lik_data_from_null = function(delta, null_fit, worker_state){
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

This means the report’s difficulty labels should be read as:

- “easy enough to support bootstrap-LR p-values first”
- with CI inversion being the same method but at significantly higher cost

## Main Architectural Gap

The current bootstrap worker abstraction is built around **resampling observed
units**:

- `bootstrap_sample_indices(...)`
- `load_bootstrap_sample_into_worker(...)`
- `bootstrap_subset_inference(...)`

Bootstrap-calibrated LR needs a second branch of worker logic:

- duplicate the inference object / design structure once
- keep `X`, `w`, and any matching/blocking structure fixed
- replace the response under a null-data simulator
- recompute unrestricted and constrained likelihood fits

So the biggest missing abstraction is not LR refitting. It is **null-data
simulation under each likelihood family**.

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
| `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm_cpp` | **Difficult** | Hybrid conditional-logit plus GLMM null simulation is too bespoke for a clean first rollout. |

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
| `InferenceCountKKHurdlePoissonOneLik` | `fast_hurdle_poisson_glmm_cpp` | **Difficult** | Hurdle plus GLMM structure is too custom for a broad first implementation. |
| `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined_with_var_cpp` | **Difficult** | Custom combined likelihood with unclear generic null-simulation semantics. |

### Continuous

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceContinKKOLSOneLik` | `fast_ols_with_var_cpp` | **Easy** | Gaussian null simulation is straightforward if the constrained fit exposes residual variance. |
| `InferenceContinKKGLMM` | `fast_gaussian_lmm_cpp` | **Borderline** | Parametric LMM simulation is feasible, but variance components and random effects add work. |

### Ordinal

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceOrdinalPropOddsRegr` | `fast_ordinal_regression_cpp` | **Borderline** | Generative simulation is possible via category probabilities, but requires family-specific ordinal sampling code. |
| `InferenceOrdinalOrderedProbitRegr` | `fast_ordinal_probit_regression_cpp` | **Borderline** | Same issue with probit thresholds. |
| `InferenceOrdinalCauchitRegr` | `fast_ordinal_cauchit_regression_cpp` | **Borderline** | Same as above. |
| `InferenceOrdinalCloglogRegr` | `fast_ordinal_cloglog_regression_cpp` | **Borderline** | Same as above. |
| `InferenceOrdinalKKGLMM` | `fast_ordinal_glmm_cpp` | **Difficult** | Ordinal thresholds plus random effects move this out of the easy bootstrap-LR tier. |

### Proportion

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferencePropBetaRegr` | `fast_beta_regression_cpp` | **Borderline** | Beta simulation is feasible, but the mean/precision parameterization needs family-specific handling. |
| `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta_cpp` | **Borderline** | Still generative, but three-part mixture simulation is more bespoke. |
| `InferencePropKKGLMM` | `fast_logistic_glmm_cpp` | **Borderline** | Same random-effects simulation issue as the other GLMM families. |

### Survival

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceSurvivalCoxPHRegr` | `fast_coxph_regression_cpp` | **Difficult** | The likelihood path is partial likelihood, not a full generative survival model. Null simulation requires extra baseline-hazard and censoring choices. |
| `InferenceSurvivalStratCoxPHRegr` | `fast_stratified_coxph_regression_cpp` | **Difficult** | Same Cox issue with added strata-specific baseline structure. |
| `InferenceSurvivalKKLWACoxOneLik` | `fast_coxph_regression_cpp` | **Difficult** | Combined-design Cox path still lacks a natural generic full-data simulator. |
| `InferenceSurvivalKKStratCoxOneLik` | `fast_stratified_coxph_regression_cpp` | **Difficult** | Same issue with stratification. |
| `InferenceSurvivalWeibullRegr` | `fast_weibull_regression_cpp` | **Borderline** | Full parametric survival simulation is possible, but event-time and censoring generation add family-specific work. |
| `InferenceSurvivalDepCensTransformRegr` | `fast_dep_cens_transform_optim_cpp` | **Difficult** | Highly bespoke joint event/censoring structure. |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty_cpp` | **Difficult** | Frailty simulation plus censoring plus repeated null refits is too heavy for a first generic rollout. |
| `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_clayton_weibull_aft_optim_cpp` | **Difficult** | Copula-based survival simulation is possible in principle but too bespoke package-wide. |

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
- the current LR plumbing in `InferenceAsympLik` can be reused almost directly

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

## Difficult Tier

These are the weakest candidates for a common bootstrap-LR rollout:

- Cox / stratified Cox / combined-design Cox
- the dependent-censoring transform path
- frailty survival models
- copula survival models
- hurdle-GLMM and other custom hybrid combined-likelihood paths

The recurring issue is that the existing likelihood test does not come with a
simple generic full-data simulator under the null.

In particular, the Cox paths are difficult because partial likelihood is not by
itself a complete generative model for survival times and censoring.

## Recommended Implementation Plan

### Phase 1: Base-Class API

Add a generic interface in `InferenceAsympLik`:

1. `compute_lik_ratio_bootstrap_two_sided_pval(delta = 0, B = 199, ...)`
2. optionally `compute_lik_ratio_bootstrap_confidence_interval(...)`
3. subclass hooks like:
   - `supports_lik_ratio_bootstrap()`
   - `simulate_under_lik_null(...)`
   - possibly `create_lik_ratio_bootstrap_worker_state(...)`

This is a moderate framework change, but it is localized and conceptually
clean.

### Phase 2: Easy Generative Families

Implement first for:

1. binary logit
2. binary probit
3. Poisson
4. modified-Poisson incidence wrappers
5. log-binomial
6. identity-binomial
7. Gaussian OLS one-likelihood path

This gives a large fraction of the package’s most standard LR paths with
relatively little model-specific simulation complexity.

### Phase 3: Selected Borderline Families

Only after the base framework is working:

1. negative binomial
2. ordinal fixed-effects models
3. Weibull
4. selected GLMMs
5. maybe zero-inflated / hurdle paths if they are worth the runtime

I would defer Cox and the more exotic survival hybrids.

## Bottom Line

Bootstrap-calibrated LR testing is a realistic extension for this repo, but the
reason is **not** that the current bootstrap is already enough. The existing
bootstrap is mainly observed-data resampling, while LR calibration needs
**parametric null simulation**.

So the real summary is:

- **Easy**: simple generative GLM-like paths and Gaussian OLS
- **Borderline**: ordinal, negative binomial, beta, Weibull, GLMM, and some
  mixture families
- **Difficult**: Cox, frailty, copula, dependent-censoring, and custom hybrid
  combined-likelihood paths

That makes bootstrap-calibrated LR a good candidate for a **selective rollout**
rather than a simultaneous all-path implementation.
