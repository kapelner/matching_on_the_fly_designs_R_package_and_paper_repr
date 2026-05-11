# Median-Bias-Correction Report For Likelihood Paths

## Scope

This report covers the package’s **likelihood-backed inference paths**, i.e. the
concrete classes participating in the `get_likelihood_test_spec()` surface used
for likelihood-ratio, score, gradient, and Wald-style inference.

The question is:

> How hard would it be to add **median bias correction** to the package’s
> likelihood-based inference paths?

Here “median bias correction” means the broad class of higher-order
bias-reduction methods that aim to make the estimator approximately
median-unbiased rather than mean-unbiased. In adjusted-score form, the target is
typically an estimator `theta_tilde` satisfying

```text
P(theta_tilde <= theta) = 1/2 + smaller-order error.
```

This is different from Firth / mean bias reduction, which targets the mean bias
of the full estimator.

I use three engineering labels:

- **Easy**: a practical median-bias-correction implementation is a realistic
  target.
- **Borderline**: conceptually plausible, but would require substantial
  model-specific higher-order derivation.
- **Difficult**: too mixture-, quadrature-, frailty-, copula-, or
  custom-likelihood-heavy for a practical common implementation.

This is an implementation audit, not a statement about theoretical existence in
all regular models.

## Short Answer

Median bias correction is usually **harder than Firth** for a general package
implementation.

Why:

- Firth often comes with a clean mean-bias-reduction / Jeffreys-penalty
  interpretation
- median bias correction usually does **not** have as simple a universal
  penalized-likelihood form
- the adjusted-score terms are often more model-specific

So while the goal is attractive, especially for skewed finite-sample
distributions, the implementation burden is typically higher than mean bias
reduction for the same family.

For this repo, median bias correction is most plausible in the same
smooth fixed-effects likelihood families that are already plausible for Firth,
but usually one notch harder.

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

The ordinary MLE solves

```text
U(theta) = 0.
```

Mean bias reduction changes this to

```text
U(theta) + A_mean(theta) = 0
```

to remove the leading mean bias term.

Median bias reduction instead uses a different adjusted score

```text
U(theta) + A_med(theta) = 0
```

where `A_med(theta)` is chosen so that the resulting estimator is better
centered in the **median** sense rather than the **mean** sense.

Operationally:

- mean bias reduction targets `E(theta_hat - theta)`
- median bias reduction targets the asymmetry of the estimator’s sampling
  distribution around the true value

So median bias correction is especially relevant when finite-sample skewness is
important.

## What Can Be Centralized

The package can centralize the outer optimization / caching / testing plumbing
in the same way as for Firth-style work:

1. define a family-specific adjusted score
2. solve the adjusted estimating equation
3. cache the corrected fit
4. decide what testing / CI machinery is supported

But unlike ridge or many Firth implementations, there is no obvious generic
package-wide penalty to plug into the current likelihood optimization stack.

So the reusable part is mostly:

- fit orchestration
- caching
- warm starts
- exposure of corrected estimates and maybe corrected Wald inference

The actual median-bias adjustment remains family-specific.

## What Is Missing

Each supported family would need either:

- an adjusted-score derivation for median bias reduction, or
- a documented equivalent correction formula for the target parameters

This means the hard part is:

- higher-order derivatives / expansions
- parameter-block-specific correction terms
- deciding whether the correction is for the full parameter or only the
  treatment component

That is a steeper ask than ordinary mean bias reduction in many families.

## Effect On Testing And Confidence Intervals

This is similar in spirit to the explicit bias-correction report, but with a
slightly different emphasis.

### Wald

If the package computes a median-bias-corrected estimate and an accompanying
variance estimate, then corrected Wald inference is the easiest path:

```text
Z_med = (theta_tilde - delta) / se_tilde.
```

So:

- **Wald p-values** can be updated most naturally
- **Wald confidence intervals** can be updated most naturally

This is the cleanest first inferential companion.

### Score

If the method is implemented as an adjusted-score estimator, then in principle a
median-bias-corrected score test might be definable for some families. But there
is no generic package-wide shortcut.

So a first implementation should assume:

- ordinary score tests remain unchanged unless a family-specific corrected score
  theory is added

### Likelihood Ratio

There is generally no simple universal median-bias-corrected likelihood-ratio
statistic analogous to “just use the same LR with a corrected point estimate”.

So unless a family-specific likelihood modification is known:

- ordinary LR tests remain unchanged
- ordinary LR-inverted CIs remain unchanged

### Gradient

The same logic applies:

- no generic package-wide median-bias-corrected gradient test
- ordinary gradient p-values / CIs would usually remain unchanged

### Practical Consequence

The most realistic first implementation is:

- corrected estimate
- corrected Wald inference
- unchanged score / LR / gradient methods

That avoids claiming a level of inferential coherence that median bias
correction alone does not automatically provide.

## Audit Table

### Incidence

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceIncidLogRegr` | `fast_logistic_regression_cpp` | **Borderline** | Standard binary GLM and a plausible target, but median-bias adjustments are more bespoke than Firth. |
| `InferenceIncidProbitRegr` | `fast_ordinal_probit_regression_cpp` in the 2-category case | **Borderline** | Smooth binary probit likelihood, but no clean generic plug-in correction. |
| `InferenceIncidModifiedPoisson` | `fast_poisson_regression_cpp` | **Borderline** | Same Poisson companion-likelihood structure as ordinary Poisson, but median bias is not as standard as mean bias here. |
| `InferenceIncidKKModifiedPoisson` | `fast_poisson_regression_cpp` | **Borderline** | Same as above. |
| `InferenceIncidLogBinomial` | `fast_log_binomial_regression_cpp` | **Borderline** | Plausible but delicate because of boundary-sensitive geometry. |
| `InferenceIncidBinomialIdentityRiskDiff` | `fast_identity_binomial_regression_cpp` | **Borderline** | Same issue with an awkward parameterization. |
| `InferenceIncidKKClogitOneLik` | stacked conditional-logistic + reservoir-logistic path via `fast_logistic_regression_with_var_cpp` | **Difficult** | Custom combined likelihood with no off-the-shelf median-bias-reduction recipe. |
| `InferenceIncidKKGLMM` | `fast_logistic_glmm_cpp` | **Difficult** | GLMM plus quadrature plus variance components is too bespoke for an early rollout. |
| `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm_cpp` | **Difficult** | Hybrid conditional-logit plus GLMM structure is too specialized. |

### Count

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceCountPoisson` | `fast_poisson_regression_cpp` | **Borderline** | Smooth and standard, but median-bias-reduction formulas are less canonical than for mean bias reduction. |
| `InferenceCountRobustPoisson` | `fast_poisson_regression_cpp` | **Borderline** | Same likelihood path as Poisson. |
| `InferenceCountQuasiPoisson` | `fast_poisson_regression_cpp` | **Borderline** | Same comment as above. |
| `InferenceCountNegBin` | `fast_neg_bin_cpp`, `fast_neg_bin_with_var_cpp` | **Difficult** | Dispersion parameter plus higher-order asymmetry makes this materially bespoke. |
| `InferenceCountZeroInflatedPoisson` | `fast_zero_augmented_poisson_cpp` | **Difficult** | Mixture structure makes median-bias derivation impractical as a common feature. |
| `InferenceCountHurdlePoisson` | `fast_zero_augmented_poisson_cpp` | **Difficult** | Same issue with hurdle structure. |
| `InferenceCountZeroInflatedNegBin` | `fast_zinb_cpp` | **Difficult** | Mixture plus dispersion is too specialized. |
| `InferenceCountHurdleNegBin` | `fast_hurdle_negbin_cpp` | **Difficult** | Same issue as above. |
| `InferenceCountKKGLMM` | `fast_poisson_glmm_cpp` | **Difficult** | Random effects and quadrature make this too custom for an early median-bias rollout. |
| `InferenceCountKKHurdlePoissonOneLik` | `fast_hurdle_poisson_glmm_cpp` | **Difficult** | Hurdle plus GLMM is too bespoke. |
| `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined_with_var_cpp` | **Difficult** | Custom combined likelihood with nonstandard nuisance geometry. |

### Continuous

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceContinKKOLSOneLik` | `fast_ols_with_var_cpp` | **Easy** | Gaussian linear-model structure is the cleanest target. |
| `InferenceContinKKGLMM` | `fast_gaussian_lmm_cpp` | **Difficult** | Variance components and random effects make median-bias derivations substantially harder. |

### Ordinal

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceOrdinalPropOddsRegr` | `fast_ordinal_regression_cpp` | **Difficult** | Threshold parameters and higher-order link-specific asymmetry make this bespoke. |
| `InferenceOrdinalOrderedProbitRegr` | `fast_ordinal_probit_regression_cpp` | **Difficult** | Same issue with a probit link. |
| `InferenceOrdinalCauchitRegr` | `fast_ordinal_cauchit_regression_cpp` | **Difficult** | Same cutpoint problem plus link-specific higher-order terms. |
| `InferenceOrdinalCloglogRegr` | `fast_ordinal_cloglog_regression_cpp` | **Difficult** | Same issue as above. |
| `InferenceOrdinalKKGLMM` | `fast_ordinal_glmm_cpp` | **Difficult** | Thresholds plus random effects plus quadrature are too specialized. |

### Proportion

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferencePropBetaRegr` | `fast_beta_regression_cpp` | **Difficult** | Mean-plus-precision structure makes median-bias correction materially more bespoke than mean-bias correction. |
| `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta_cpp` | **Difficult** | Three-part mixture is too specialized. |
| `InferencePropKKGLMM` | `fast_logistic_glmm_cpp` | **Difficult** | Same quadrature GLMM issue as the incidence counterpart. |

### Survival

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceSurvivalCoxPHRegr` | `fast_coxph_regression_cpp` | **Difficult** | Partial-likelihood median-bias reduction is not a practical common target here. |
| `InferenceSurvivalStratCoxPHRegr` | `fast_stratified_coxph_regression_cpp` | **Difficult** | Same issue with stratification. |
| `InferenceSurvivalKKLWACoxOneLik` | `fast_coxph_regression_cpp` | **Difficult** | Combined-design Cox path is too bespoke. |
| `InferenceSurvivalKKStratCoxOneLik` | `fast_stratified_coxph_regression_cpp` | **Difficult** | Same with stratification. |
| `InferenceSurvivalWeibullRegr` | `fast_weibull_regression_cpp` | **Borderline** | Smooth parametric survival model, but substantially more bespoke than Gaussian and basic GLMs. |
| `InferenceSurvivalDepCensTransformRegr` | `fast_dep_cens_transform_optim_cpp` | **Difficult** | Coupled event/censoring structure is too custom. |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty_cpp` | **Difficult** | Frailty integration and variance terms are too bespoke. |
| `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_clayton_weibull_aft_optim_cpp` | **Difficult** | Copula dependence plus survival margins is too specialized. |

## Easy Tier

These are the only clearly realistic first targets:

1. `InferenceContinKKOLSOneLik`

That is the cleanest case because Gaussian likelihood geometry is simplest and
the distinction between mean and median centering is least burdensome to
operationalize.

## Borderline Tier

These are the best plausible non-Gaussian targets:

- ordinary logit
- ordinary probit
- Poisson and companion Poisson paths
- log-binomial
- identity-binomial
- Weibull

These are the families where a serious family-specific effort could plausibly
pay off, but I would not call them easy.

## Difficult Tier

These are poor targets for a common implementation:

- negative binomial
- ordinal fixed-effects models
- beta regression
- all mixture likelihoods
- all GLMM paths
- Cox / stratified Cox
- frailty, copula, and dependent-censoring survival paths
- custom combined likelihood hybrids

For these families, median-bias correction would quickly become a model-specific
research project.

## Recommended Implementation Plan

### Phase 1: Decide Whether Median Bias Is Worth Pursuing At All

Because this is usually harder than Firth, the package should first decide
whether the inferential gain is worth the additional family-specific work.

If the real goal is:

- mean bias reduction
- separation resistance
- simple smooth corrected likelihoods

then Firth is usually the better first project.

### Phase 2: Expose Explicit Median-Bias-Corrected Wald Support

If median bias correction is pursued, the cleanest first API is:

1. `compute_median_bias_corrected_estimate()`
2. `compute_median_bias_corrected_wald_two_sided_pval()`
3. `compute_median_bias_corrected_wald_confidence_interval()`

That avoids overstating the extent to which score / LR / gradient inference has
been corrected.

### Phase 3: Start With The Simplest Smooth Families

If implemented at all, start with:

1. Gaussian OLS one-likelihood path
2. perhaps ordinary logit
3. perhaps ordinary probit
4. perhaps Poisson

Only after that should the package consider the harder families.

## Bottom Line

Median bias correction is conceptually appealing, but for this codebase it is
generally **harder than Firth** and much less likely to admit a clean common
implementation across the likelihood-backed inference paths.

The most realistic interpretation is:

- family-specific corrected point estimates
- corrected Wald inference where available
- unchanged score / LR / gradient methods unless separately developed

So if the package wants a broadly useful higher-order likelihood improvement,
median bias correction is probably a **later-stage selective extension**, not an
early package-wide project.
