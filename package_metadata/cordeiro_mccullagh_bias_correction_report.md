# Cordeiro-McCullagh Bias-Correction Report For Likelihood Paths

## Scope

This report covers the package’s **likelihood-backed inference paths**, i.e. the
concrete classes participating in the `get_likelihood_test_spec()` surface used
for likelihood-ratio, score, gradient, and Wald-style inference.

The question is:

> How hard would it be to add **Cordeiro-McCullagh type explicit bias
> correction** to the package’s likelihood-based inference paths?

Here the target is the usual first-order analytic bias correction to the MLE:

```text
theta_tilde = theta_hat - b_hat / n
```

where:

- `theta_hat` is the ordinary MLE
- `b_hat / n` is an estimated `O(n^{-1})` bias term
- `theta_tilde` is the bias-corrected estimator

This is **not** a penalized likelihood and **not** an adjusted-score method. The
ordinary likelihood fit is computed first, and then the estimate is shifted by
an estimated bias term.

I use three engineering labels:

- **Easy**: explicit analytic bias correction is a realistic implementation
  target.
- **Borderline**: conceptually plausible, but would require substantial
  model-specific higher-order derivations.
- **Difficult**: too mixture-, quadrature-, frailty-, copula-, or
  custom-likelihood-heavy for a practical common implementation.

This is an implementation audit, not a claim about abstract mathematical
existence.

## Short Answer

This approach is attractive for the package because it is often **lighter-touch
than Firth**:

- it keeps the ordinary MLE optimization unchanged
- it does not require a penalized objective
- it does not require changing `fit_null(...)` optimizers just to get a
  corrected point estimate

But it is also weaker and messier from an inference perspective:

- it naturally gives a **bias-corrected point estimate**
- it does **not automatically define** a new coherent score/LR/gradient test
- it does **not automatically define** a new confidence interval method
- it is most naturally paired with **Wald-style inference**

So the best use case in this repo would be:

1. compute the ordinary likelihood fit
2. compute an analytic bias estimate for the treatment coefficient
3. report a bias-corrected treatment estimate
4. initially keep score / gradient / LR tests tied to the uncorrected
   likelihood, unless a family-specific corrected-testing theory is supplied

## Core Math

Let `theta` denote the full parameter vector, with log-likelihood

```text
l(theta),
```

score

```text
U(theta) = d l(theta) / d theta,
```

and information matrix

```text
I(theta).
```

Under standard regularity conditions, the MLE has expansion

```text
E(theta_hat - theta) = b(theta) / n + O(n^{-2}),
```

for a model-specific bias function `b(theta)`.

The Cordeiro-McCullagh idea is:

1. fit the ordinary MLE `theta_hat`
2. estimate the leading bias term `b(theta_hat) / n`
3. define

```text
theta_tilde = theta_hat - b(theta_hat) / n
```

The resulting estimator typically has smaller mean bias:

```text
E(theta_tilde - theta) = O(n^{-2})
```

if the bias estimate is correct to first order.

Unlike Firth:

- the likelihood is not changed
- the score equation is not changed
- the null-constrained fit is not changed

Instead, the correction is applied **after fitting**.

## What Can Be Centralized

The package already has most of what is needed to support the **ordinary** MLE
fit:

- `compute_estimate()`
- access to fitted parameters
- `get_likelihood_test_spec()`
- access to score and information in many paths

So a generic design is plausible:

1. compute the ordinary full fit
2. ask the subclass for the bias estimate for the full parameter or treatment
   component
3. store both:
   - ordinary estimate
   - bias-corrected estimate

Conceptually, a family hook could look like:

```r
get_mle_bias_correction = function(spec, fit){
  list(
    bias_full = ...,
    bias_treatment = ...
  )
}
```

or more simply:

```r
get_treatment_bias_correction = function(spec, fit){
  ...
}
```

Then the base layer could expose:

- `compute_estimate()` as the corrected estimate for supported families, or
- a separate method like `compute_bias_corrected_estimate()`

The second option is safer, because it avoids silently changing the meaning of
existing inference objects.

## What Is Missing

The hard part is the same one that appears in other higher-order corrections:

- the leading bias term is **family-specific**
- it often requires higher-order derivatives and expectations
- it is usually easier in standard GLM-like models than in mixtures, GLMMs, or
  custom combined likelihoods

So while the package can centralize the storage and dispatch of a bias
correction, it cannot centralize the actual formula in a model-agnostic way.

## Effect On Testing And Confidence Intervals

This is the most important design issue for this method.

### Wald Paths

This method fits most naturally with Wald inference.

If the treatment estimate is corrected from

```text
beta_hat_T
```

to

```text
beta_tilde_T = beta_hat_T - bias_hat_T,
```

then the simplest Wald statistic becomes

```text
Z_tilde = (beta_tilde_T - delta) / se_tilde
```

where `se_tilde` could be handled in one of three ways:

1. **plug-in ordinary SE**
   use the same information-based SE as the uncorrected MLE
2. **bias-corrected asymptotic variance**
   derive a refined variance for the corrected estimator
3. **delta-method / sandwich refinement**
   treat the corrected estimator as a smooth function of the MLE

For a first implementation, option 1 is the most realistic.

So:

- **Wald p-values** can be updated immediately using the corrected estimate
- **Wald CIs** can also be updated immediately
- these would be the most natural inferential companion to the corrected
  estimator

### Score Paths

The ordinary score test is based on the **uncorrected** score at the null:

```text
U(theta_0)
```

with information under the null.

An explicit post-fit bias correction does **not** define a new score equation.
So unless you derive a corrected score test separately, the ordinary score test
remains the natural score test.

That means:

- the **score p-value** would usually stay tied to the original likelihood
- the **score CI** obtained by inversion would also stay tied to the original
  score test

So score inference would generally **not** align automatically with the
bias-corrected point estimate.

### Likelihood-Ratio Paths

The likelihood-ratio test depends on:

```text
LR(delta) = 2 { l(theta_hat) - l(theta_hat_delta) }.
```

Because the Cordeiro-McCullagh correction does not change the likelihood, there
is no automatic corrected LR statistic.

You could still report the ordinary LR test, but it would be a test for the
ordinary likelihood fit, not for a corrected-likelihood estimator.

So:

- **LR p-values** would ordinarily remain unchanged
- **LR CIs** obtained by inversion would ordinarily remain unchanged

If the package wants LR-based inference to reflect the correction, that becomes
a different project, closer to modified profile or Bartlett-style work.

### Gradient Paths

The gradient test in the repo is also based on the restricted score and the
ordinary estimate under the current likelihood framework.

A post-fit bias correction does not automatically define a corrected gradient
test either.

So:

- **gradient p-values** would ordinarily remain unchanged
- **gradient CIs** would ordinarily remain unchanged

### Practical Consequence

This method can easily create an inferential mismatch:

- corrected point estimate
- uncorrected score / LR / gradient tests and intervals

That is not necessarily wrong, but it needs to be explicit in the API and
documentation.

## Recommended Inference Policy

If this method is added, the cleanest first policy would be:

1. expose the corrected estimator explicitly
2. use corrected Wald inference for supported families
3. leave score / gradient / LR inference unchanged unless a family-specific
   corrected theory is implemented

So the package might expose something like:

- `compute_bias_corrected_estimate()`
- `compute_bias_corrected_wald_two_sided_pval()`
- `compute_bias_corrected_wald_confidence_interval()`

while leaving:

- `compute_score_two_sided_pval()`
- `compute_gradient_two_sided_pval()`
- `compute_lik_ratio_two_sided_pval()`

as ordinary-likelihood methods.

That is much cleaner than pretending all testing paths have been corrected when
only the point estimate has.

## Audit Table

### Incidence

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceIncidLogRegr` | `fast_logistic_regression_cpp` | **Easy** | Standard binary GLM; explicit first-order bias formulas are realistic. |
| `InferenceIncidProbitRegr` | `fast_ordinal_probit_regression_cpp` in the 2-category case | **Easy** | Smooth binary probit likelihood with standard fixed-effects structure. |
| `InferenceIncidModifiedPoisson` | `fast_poisson_regression_cpp` | **Easy** | Same Poisson companion-likelihood structure as ordinary Poisson. |
| `InferenceIncidKKModifiedPoisson` | `fast_poisson_regression_cpp` | **Easy** | Same as above. |
| `InferenceIncidLogBinomial` | `fast_log_binomial_regression_cpp` | **Borderline** | Still a regular GLM-type model, but less canonical and more delicate near boundaries. |
| `InferenceIncidBinomialIdentityRiskDiff` | `fast_identity_binomial_regression_cpp` | **Borderline** | Bias formulas are plausible but awkward because of boundary-sensitive parameterization. |
| `InferenceIncidKKClogitOneLik` | stacked conditional-logistic + reservoir-logistic path via `fast_logistic_regression_with_var_cpp` | **Borderline** | Custom combined likelihood with no off-the-shelf bias formula. |
| `InferenceIncidKKGLMM` | `fast_logistic_glmm_cpp` | **Difficult** | Quadrature GLMM likelihood with variance components is not a clean first-order bias-correction target. |
| `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm_cpp` | **Difficult** | Hybrid conditional-logit plus GLMM structure is too bespoke. |

### Count

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceCountPoisson` | `fast_poisson_regression_cpp` | **Easy** | Canonical Poisson GLM and a natural first target. |
| `InferenceCountRobustPoisson` | `fast_poisson_regression_cpp` | **Easy** | The likelihood fit is still ordinary Poisson. |
| `InferenceCountQuasiPoisson` | `fast_poisson_regression_cpp` | **Easy** | Same LR/likelihood path as Poisson. |
| `InferenceCountNegBin` | `fast_neg_bin_cpp`, `fast_neg_bin_with_var_cpp` | **Borderline** | Feasible in principle, but dispersion-parameter bias terms add substantial model-specific work. |
| `InferenceCountZeroInflatedPoisson` | `fast_zero_augmented_poisson_cpp` | **Difficult** | Mixture structure makes the explicit bias term highly bespoke. |
| `InferenceCountHurdlePoisson` | `fast_zero_augmented_poisson_cpp` | **Difficult** | Same issue with hurdle/truncation structure. |
| `InferenceCountZeroInflatedNegBin` | `fast_zinb_cpp` | **Difficult** | Mixture plus dispersion parameter is too specialized. |
| `InferenceCountHurdleNegBin` | `fast_hurdle_negbin_cpp` | **Difficult** | Same issue as above. |
| `InferenceCountKKGLMM` | `fast_poisson_glmm_cpp` | **Difficult** | GLMM plus quadrature plus variance component is not a realistic first wave. |
| `InferenceCountKKHurdlePoissonOneLik` | `fast_hurdle_poisson_glmm_cpp` | **Difficult** | Hurdle plus GLMM structure is too bespoke. |
| `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined_with_var_cpp` | **Difficult** | Custom combined likelihood with nonstandard nuisance geometry. |

### Continuous

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceContinKKOLSOneLik` | `fast_ols_with_var_cpp` | **Easy** | Gaussian linear-model setting is a natural explicit bias-correction target. |
| `InferenceContinKKGLMM` | `fast_gaussian_lmm_cpp` | **Borderline** | Possible in principle, but variance components make the first-order bias term more bespoke. |

### Ordinal

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceOrdinalPropOddsRegr` | `fast_ordinal_regression_cpp` | **Borderline** | Fixed-effects ordinal model is smooth, but threshold parameters make bias formulas heavier. |
| `InferenceOrdinalOrderedProbitRegr` | `fast_ordinal_probit_regression_cpp` | **Borderline** | Same threshold issue with a probit link. |
| `InferenceOrdinalCauchitRegr` | `fast_ordinal_cauchit_regression_cpp` | **Borderline** | Same with a different link-specific likelihood geometry. |
| `InferenceOrdinalCloglogRegr` | `fast_ordinal_cloglog_regression_cpp` | **Borderline** | Same issue as above. |
| `InferenceOrdinalKKGLMM` | `fast_ordinal_glmm_cpp` | **Difficult** | Thresholds plus random effects plus quadrature are too bespoke. |

### Proportion

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferencePropBetaRegr` | `fast_beta_regression_cpp` | **Borderline** | Smooth model, but mean-plus-precision bias formulas are model-specific and more involved than basic GLMs. |
| `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta_cpp` | **Difficult** | Three-part mixture structure is too specialized. |
| `InferencePropKKGLMM` | `fast_logistic_glmm_cpp` | **Difficult** | Same quadrature GLMM issue as the incidence counterpart. |

### Survival

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceSurvivalCoxPHRegr` | `fast_coxph_regression_cpp` | **Borderline** | Cox partial likelihood bias correction is plausible, but not as straightforward as standard GLMs. |
| `InferenceSurvivalStratCoxPHRegr` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same issue with added stratification structure. |
| `InferenceSurvivalKKLWACoxOneLik` | `fast_coxph_regression_cpp` | **Borderline** | Combined-design Cox path is plausible but bespoke. |
| `InferenceSurvivalKKStratCoxOneLik` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same with stratification. |
| `InferenceSurvivalWeibullRegr` | `fast_weibull_regression_cpp` | **Borderline** | Smooth parametric survival model, but bias formulas are materially more bespoke than basic GLMs. |
| `InferenceSurvivalDepCensTransformRegr` | `fast_dep_cens_transform_optim_cpp` | **Difficult** | Coupled event/censoring structure is too custom. |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty_cpp` | **Difficult** | Frailty integration and variance terms make this a poor target for simple explicit bias correction. |
| `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_clayton_weibull_aft_optim_cpp` | **Difficult** | Copula dependence plus survival margins is too bespoke. |

## Easy Tier

These are the most realistic first targets:

1. `InferenceIncidLogRegr`
2. `InferenceIncidProbitRegr`
3. `InferenceCountPoisson`
4. `InferenceCountRobustPoisson`
5. `InferenceCountQuasiPoisson`
6. `InferenceIncidModifiedPoisson`
7. `InferenceIncidKKModifiedPoisson`
8. `InferenceContinKKOLSOneLik`

These are best because:

- the MLE fit already exists and is stable
- first-order bias expansions are most plausible here
- corrected Wald inference is easy to layer on top

## Borderline Tier

These are plausible but would require dedicated formulas:

- log-binomial
- identity-binomial
- negative binomial
- ordinal fixed-effects models
- beta regression
- Cox / stratified Cox paths
- Weibull
- Gaussian LMM
- combined conditional-logistic path

These are reasonable second-wave candidates if the package wants explicit
bias-corrected point estimates beyond the easiest GLM families.

## Difficult Tier

These are poor targets for a common implementation:

- zero-inflated and hurdle count models
- zero/one-inflated beta
- all quadrature GLMM paths except perhaps as a separate research project
- frailty models
- copula models
- dependent-censoring transform likelihood
- custom hybrid combined-likelihood paths

For these families, the explicit bias term quickly becomes too bespoke to treat
as a common package feature.

## Recommended Implementation Plan

### Phase 1: Explicit Estimate API

Do **not** silently redefine the existing estimate first.

Instead, add explicit methods such as:

1. `compute_bias_corrected_estimate()`
2. `compute_bias_corrected_wald_two_sided_pval()`
3. `compute_bias_corrected_wald_confidence_interval()`

This preserves the current meaning of ordinary MLE-based methods and avoids
mixing corrected and uncorrected inferential semantics by accident.

### Phase 2: Easy Families First

Implement first for:

1. ordinary logit
2. ordinary probit
3. Poisson
4. modified-Poisson incidence wrappers
5. Gaussian OLS one-likelihood path

Optionally extend next to:

6. log-binomial
7. identity-binomial

### Phase 3: Decide On Inference Policy

The cleanest initial policy is:

- corrected estimate
- corrected Wald inference
- ordinary score / gradient / LR inference left unchanged

If that mismatch is judged undesirable, then this approach may not be worth
using for families where the package relies heavily on score/LR/gradient-based
inference.

### Phase 4: Selective Expansion

Only after that foundation is stable, consider:

1. negative binomial
2. ordinal fixed-effects models
3. beta regression
4. Weibull
5. perhaps selected Cox paths

I would not plan on a package-wide rollout to the mixture, GLMM, frailty,
copula, and custom-hybrid paths.

## Bottom Line

Cordeiro-McCullagh-type explicit bias correction is a realistic way to improve
point estimation in the simpler likelihood paths because it leaves the ordinary
likelihood optimization untouched.

Its main limitation is inferential coherence:

- it naturally improves the **point estimate**
- it can be paired cleanly with **Wald inference**
- it does **not** automatically redefine score, gradient, or likelihood-ratio
  tests and inverted CIs

So for this repo, the method is best viewed as:

- a plausible extension for simple GLM-like and Gaussian likelihood paths
- a less natural fit for the package’s score/LR/gradient-heavy inference logic
- a selective feature, not a simultaneous all-path implementation
