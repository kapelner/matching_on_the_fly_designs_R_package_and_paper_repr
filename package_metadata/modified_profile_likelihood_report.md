# Modified Profile Likelihood / Cox-Reid Report For Likelihood Paths

## Scope

This report covers the package’s **likelihood-backed inference paths**, i.e. the
concrete classes participating in the `get_likelihood_test_spec()` surface used
for likelihood-ratio p-values and inverted confidence intervals.

The question is:

> How hard would it be to add **modified profile likelihood** inference for the
> treatment effect across the package’s likelihood-based inference paths?

Here “modified profile likelihood” includes the Cox-Reid style family of
adjustments where:

- the treatment effect is the scalar parameter of interest
- the remaining parameters are nuisance parameters
- inference is based on an adjusted profile likelihood rather than the ordinary
  profile likelihood

I use three engineering labels:

- **Easy**: a modified-profile implementation is a natural fit for the path.
- **Borderline**: conceptually plausible, but would require substantial
  family-specific derivation or plumbing.
- **Difficult**: too mixture-, quadrature-, frailty-, copula-, or
  custom-likelihood-heavy for a practical common implementation.

This is an implementation audit, not a claim about abstract mathematical
existence.

## Short Answer

Modified profile likelihood is most attractive when the package wants to improve
**inference for the treatment effect** rather than globally redefining the full
estimator.

Compared with the other small-sample corrections:

- unlike **Firth**, it does not primarily target full-parameter bias reduction
- unlike **Bartlett**, it does not just rescale the LR statistic afterward
- unlike **bootstrap-calibrated LR**, it does not replace asymptotics with
  simulation

Instead, it changes the likelihood used to profile the treatment effect by
adjusting for nuisance-parameter curvature.

That makes it especially appealing for paths with meaningful nuisance structure:

- dispersion parameters
- precision parameters
- threshold / cutpoint parameters
- variance components
- parametric survival nuisance structure

## Core Math

Split the parameter into:

- `psi`: the treatment effect
- `lambda`: nuisance parameters

with full log-likelihood

```text
l(psi, lambda)
```

The ordinary profile log-likelihood is

```text
l_p(psi) = l(psi, lambda_hat_psi)
```

where `lambda_hat_psi` is the nuisance MLE at fixed `psi`.

The ordinary likelihood-ratio statistic for testing `psi = psi_0` is

```text
LR(psi_0) = 2 { l_p(psi_hat) - l_p(psi_0) }.
```

The Cox-Reid style modified profile likelihood replaces `l_p` by

```text
l_MP(psi) = l(psi, lambda_hat_psi) - 1/2 log | I_lambda lambda(psi, lambda_hat_psi) |
```

up to conventions about expected vs observed nuisance information and additive
constants that do not depend on `psi`.

Here:

```text
I_lambda lambda(psi, lambda)
```

is the nuisance-information block.

The corresponding modified-profile likelihood-ratio statistic is

```text
LR_MP(psi_0) = 2 { l_MP(psi_hat_MP) - l_MP(psi_0) }.
```

In practice, for this repo, the key operational point is:

- the package already knows how to refit the model at fixed treatment effect
- what is missing is access to the nuisance-information determinant along that
  constrained path

## What Can Be Centralized

The generic outer plumbing can live in `InferenceAsympLik`.

The existing likelihood-test framework already has:

- `get_likelihood_test_spec()`
- constrained null refits at fixed `delta`
- cached full and null negative log-likelihoods
- likelihood-ratio p-value computation
- LR-based CI inversion

A modified-profile extension could therefore reuse most of that structure:

1. fit the constrained model at `delta`
2. compute the ordinary constrained log-likelihood
3. compute the nuisance-information adjustment
4. form the modified profile log-likelihood
5. use the modified LR statistic for p-values and CI inversion

So the shared part is feasible.

## What Is Missing

The hard part is family-specific access to the nuisance block.

For a clean implementation, each supported family would need a hook like:

```r
get_modified_profile_components = function(spec, fit, delta){
  list(
    profile_loglik = ...,
    nuisance_information = ...
  )
}
```

or equivalently:

```r
get_modified_profile_loglik = function(spec, fit, delta){
  ...
}
```

This means the package must know, for each family:

- which parameters are nuisance parameters
- how to extract the nuisance-information block
- whether expected or observed nuisance information is appropriate
- how to evaluate the determinant stably

That is why modified profile likelihood is naturally selective rather than
package-universal.

## Why The CI Story Is Straightforward Once The Statistic Exists

As with Bartlett and ordinary LR paths, the CI is obtained by inversion:

- ordinary LR CI:
  invert `LR(delta)`
- modified-profile CI:
  invert `LR_MP(delta)`

So once the modified-profile test statistic exists as a function of `delta`, the
confidence interval can reuse the same inversion machinery conceptually.

The extra runtime comes from evaluating the nuisance-information adjustment at
each constrained fit visited during inversion.

## Audit Table

### Incidence

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceIncidLogRegr` | `fast_logistic_regression_cpp` | **Borderline** | Possible, but this is a regular fixed-effects binary GLM where Firth or Bartlett are cleaner first improvements. |
| `InferenceIncidProbitRegr` | `fast_ordinal_probit_regression_cpp` in the 2-category case | **Borderline** | Same issue as logit: feasible, but not the highest-value first modified-profile target. |
| `InferenceIncidModifiedPoisson` | `fast_poisson_regression_cpp` | **Borderline** | Smooth and regular, but nuisance structure is modest. |
| `InferenceIncidKKModifiedPoisson` | `fast_poisson_regression_cpp` | **Borderline** | Same Poisson companion path with limited nuisance complexity. |
| `InferenceIncidLogBinomial` | `fast_log_binomial_regression_cpp` | **Borderline** | Plausible, but still a relatively regular fixed-effects GLM. |
| `InferenceIncidBinomialIdentityRiskDiff` | `fast_identity_binomial_regression_cpp` | **Borderline** | Same as above, though somewhat more awkward algebraically. |
| `InferenceIncidKKClogitOneLik` | stacked conditional-logistic + reservoir-logistic path via `fast_logistic_regression_with_var_cpp` | **Borderline** | Custom combined likelihood with meaningful nuisance structure, but not plug-and-play. |
| `InferenceIncidKKGLMM` | `fast_logistic_glmm_cpp` | **Easy** | Variance components are classic nuisance parameters for modified-profile treatment-effect inference. |
| `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm_cpp` | **Difficult** | Hybrid conditional-logit plus GLMM structure is too bespoke for a clean shared implementation. |

### Count

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceCountPoisson` | `fast_poisson_regression_cpp` | **Borderline** | Regular GLM with limited nuisance complexity. |
| `InferenceCountRobustPoisson` | `fast_poisson_regression_cpp` | **Borderline** | Same LR path as ordinary Poisson. |
| `InferenceCountQuasiPoisson` | `fast_poisson_regression_cpp` | **Borderline** | Same comment as above. |
| `InferenceCountNegBin` | `fast_neg_bin_cpp`, `fast_neg_bin_with_var_cpp` | **Easy** | Nuisance dispersion makes this one of the strongest modified-profile targets. |
| `InferenceCountZeroInflatedPoisson` | `fast_zero_augmented_poisson_cpp` | **Difficult** | Mixture structure makes clean nuisance-adjusted profiling difficult. |
| `InferenceCountHurdlePoisson` | `fast_zero_augmented_poisson_cpp` | **Difficult** | Two-part truncation structure is too bespoke. |
| `InferenceCountZeroInflatedNegBin` | `fast_zinb_cpp` | **Difficult** | Mixture plus dispersion is not a realistic common target. |
| `InferenceCountHurdleNegBin` | `fast_hurdle_negbin_cpp` | **Difficult** | Same issue with hurdle structure. |
| `InferenceCountKKGLMM` | `fast_poisson_glmm_cpp` | **Easy** | Variance-component nuisance structure makes this a natural conceptual fit. |
| `InferenceCountKKHurdlePoissonOneLik` | `fast_hurdle_poisson_glmm_cpp` | **Difficult** | Hurdle plus GLMM structure is too specialized. |
| `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined_with_var_cpp` | **Difficult** | Custom combined likelihood with nonstandard nuisance geometry. |

### Continuous

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceContinKKOLSOneLik` | `fast_ols_with_var_cpp` | **Borderline** | In principle possible, but OLS is too regular for modified profiling to be an especially valuable first gain. |
| `InferenceContinKKGLMM` | `fast_gaussian_lmm_cpp` | **Easy** | Variance components are classic nuisance parameters for this type of adjustment. |

### Ordinal

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceOrdinalPropOddsRegr` | `fast_ordinal_regression_cpp` | **Easy** | Treatment-effect profiling against threshold nuisance parameters is a strong conceptual fit. |
| `InferenceOrdinalOrderedProbitRegr` | `fast_ordinal_probit_regression_cpp` | **Easy** | Same threshold-nuisance logic with a probit link. |
| `InferenceOrdinalCauchitRegr` | `fast_ordinal_cauchit_regression_cpp` | **Easy** | Same cutpoint-based nuisance structure. |
| `InferenceOrdinalCloglogRegr` | `fast_ordinal_cloglog_regression_cpp` | **Easy** | Same as above. |
| `InferenceOrdinalKKGLMM` | `fast_ordinal_glmm_cpp` | **Easy** | Thresholds plus variance components make nuisance-adjusted profiling conceptually attractive. |

### Proportion

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferencePropBetaRegr` | `fast_beta_regression_cpp` | **Easy** | Mean-plus-precision structure makes this one of the clearest modified-profile targets. |
| `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta_cpp` | **Difficult** | Three-part mixture structure is too bespoke. |
| `InferencePropKKGLMM` | `fast_logistic_glmm_cpp` | **Easy** | Same variance-component nuisance logic as the other GLMM families. |

### Survival

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceSurvivalCoxPHRegr` | `fast_coxph_regression_cpp` | **Borderline** | Cox-Reid ideas are conceptually close, but partial likelihood makes implementation materially more bespoke. |
| `InferenceSurvivalStratCoxPHRegr` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same issue with added stratification structure. |
| `InferenceSurvivalKKLWACoxOneLik` | `fast_coxph_regression_cpp` | **Borderline** | Combined-design Cox partial likelihood is plausible but not plug-and-play. |
| `InferenceSurvivalKKStratCoxOneLik` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same with stratified risk sets. |
| `InferenceSurvivalWeibullRegr` | `fast_weibull_regression_cpp` | **Easy** | Smooth parametric survival likelihood with meaningful nuisance structure. |
| `InferenceSurvivalDepCensTransformRegr` | `fast_dep_cens_transform_optim_cpp` | **Difficult** | Coupled event/censoring parameter blocks make clean profiling very bespoke. |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty_cpp` | **Difficult** | Frailty integration moves this out of the easy modified-profile regime. |
| `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_clayton_weibull_aft_optim_cpp` | **Difficult** | Copula dependence plus parametric margins is too bespoke for a shared rollout. |

## Easy Tier

These are the strongest conceptual targets:

1. `InferenceCountNegBin`
2. `InferencePropBetaRegr`
3. `InferenceSurvivalWeibullRegr`
4. `InferenceContinKKGLMM`
5. `InferenceIncidKKGLMM`
6. `InferenceCountKKGLMM`
7. `InferenceOrdinalKKGLMM`
8. `InferencePropKKGLMM`
9. `InferenceOrdinalPropOddsRegr`
10. `InferenceOrdinalOrderedProbitRegr`
11. `InferenceOrdinalCauchitRegr`
12. `InferenceOrdinalCloglogRegr`

These are the families where nuisance-parameter adjustment is most naturally
aligned with the model structure:

- dispersion
- precision
- thresholds
- variance components
- parametric survival nuisance structure

## Borderline Tier

These can benefit, but I would not put them first:

- ordinary fixed-effects incidence GLMs
- ordinary Poisson and its robust/quasi wrappers
- OLS
- Cox / stratified Cox / combined-design Cox
- `InferenceIncidKKClogitOneLik`

The common theme is either:

- too little nuisance complexity to justify modified profiling as the first
  improvement, or
- enough bespoke structure that implementation becomes model-specific quickly

## Difficult Tier

These are poor targets for a common modified-profile rollout:

- zero-inflated and hurdle count models
- zero/one-inflated beta
- frailty models
- copula models
- dependent-censoring transform likelihood
- custom hybrids like `InferenceIncidKKClogitPlusGLMMOneLik` and
  `InferenceCountKKCPoissonOneLik`

In these families, nuisance-adjusted profiling quickly becomes a family-specific
research project rather than a reusable package feature.

## Recommended Implementation Plan

### Phase 1: Base-Class Hook

Add optional modified-profile support to `InferenceAsympLik`:

1. a new testing type, e.g. `lik_ratio_mod_profile`
2. a generic dispatcher for p-values and confidence intervals
3. a family hook such as:
   - `supports_modified_profile_likelihood()`
   - `get_modified_profile_loglik(...)`
   - or `get_modified_profile_components(...)`

This phase is a moderate refactor. The plumbing is generic, but the actual
adjustment remains family-specific.

### Phase 2: Best First Families

Implement first for:

1. negative binomial
2. beta regression
3. Weibull regression
4. one or two ordinal fixed-effects models

These give the best ratio of inferential value to implementation complexity.

### Phase 3: GLMM Families

If the approach proves valuable, extend selectively to:

1. Gaussian LMM
2. logistic GLMM
3. Poisson GLMM
4. perhaps ordinal GLMM

These are conceptually good fits, but the nuisance blocks and determinants will
take more care.

### Phase 4: Leave Exotic Families Last

I would defer:

- mixtures
- hurdles
- frailty
- copula
- dependent-censoring transforms
- custom combined-likelihood hybrids

unless there is a strong applied reason to support one specific family.

## Bottom Line

Modified profile likelihood is a good fit for this codebase when:

- the treatment effect is the parameter of interest
- nuisance parameters are substantial
- the likelihood is still smooth enough to expose a stable nuisance-information
  block

That makes the best targets:

- negative binomial
- beta regression
- Weibull regression
- ordinal fixed-effects models
- GLMM families

It is less compelling for simple fixed-effects GLMs and much less realistic for
mixtures, frailty, copula, and other custom hybrid likelihoods.

So the right strategy is a **selective rollout** with family-specific support,
not a simultaneous package-wide implementation.
