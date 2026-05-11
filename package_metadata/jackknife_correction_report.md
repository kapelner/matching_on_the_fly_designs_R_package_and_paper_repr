# Jackknife-Correction Report Across Inference Paths

## Scope

This report is **not limited to likelihood-based inference paths**. It covers
the broader package surface of inference objects that compute a treatment-effect
estimate through `compute_estimate()`.

The question is:

> How easy would it be to add **jackknife correction** to the package’s
> inference paths?

Here “jackknife correction” means the usual leave-one-out bias correction based
on the full-sample estimate `theta_hat` and the leave-one-out estimates
`theta_hat_{(-i)}`:

```text
theta_J = n * theta_hat - (n - 1) * mean(theta_hat_{(-i)})
```

equivalently

```text
theta_J = theta_hat - (n - 1) * ( mean(theta_hat_{(-i)}) - theta_hat ).
```

This is a resampling-based explicit bias correction, not a penalized likelihood
and not a modified score equation.

I use three engineering labels:

- **Easy**: a jackknife-corrected estimate is a realistic extension with the
  current architecture.
- **Borderline**: feasible, but needs special handling for design structure,
  nonsmoothness, or expensive refits.
- **Difficult**: delete-one logic is not natural for the design/model structure,
  or the computational burden / instability is too high for a clean common
  implementation.

This is an implementation audit, not a statement about abstract asymptotic
validity in every model.

## Short Answer

Jackknife correction is more broadly implementable than Firth, Bartlett,
modified profile, or bootstrap-calibrated likelihood-ratio methods because it
does **not** require:

- a likelihood
- gradients or Hessians
- constrained null refits
- parametric simulation

It only needs repeated recomputation of the treatment-effect estimate on
leave-one-out subsets.

For this repo, that matters because the package already has a major building
block in the shared bootstrap layer:

- `approximate_jackknife_distribution_beta_hat_T()`

So **estimate-level jackknife correction is partly wired already**.

The main complications are:

- what the leave-one-out unit should be in blocked, matched, clustered, or
  paired designs
- instability for nonsmooth or small-sample exact procedures
- computational cost for expensive likelihood and mixed-model paths
- inferential coherence for score / gradient / LR methods

## Core Math

If `theta_hat` is the full-sample estimate and `theta_hat_{(-i)}` is the
estimate leaving out unit `i`, define the average leave-one-out estimate

```text
theta_bar_(.) = (1/n) sum_i theta_hat_{(-i)}.
```

The classical delete-1 jackknife bias estimate is

```text
bias_hat_J = (n - 1) * (theta_bar_(.) - theta_hat).
```

The jackknife bias-corrected estimate is then

```text
theta_J = theta_hat - bias_hat_J
        = n * theta_hat - (n - 1) * theta_bar_(.).
```

This is aimed at removing leading finite-sample bias in the estimate.

Unlike Cordeiro-McCullagh or Firth:

- no analytic higher-order formula is required
- no likelihood modification is required
- the correction is entirely built from repeated refits on subsets

## Why This Is More Feasible Here Than Some Other Corrections

The package already supports resampling and subset refits through its bootstrap
infrastructure, and the shared bootstrap class already includes:

- leave-one-out estimate generation for BCa via
  `approximate_jackknife_distribution_beta_hat_T()`
- design duplication and subset inference helpers
- many inference paths with reusable worker/state patterns

So for many classes, the package is already close to being able to compute:

1. `theta_hat`
2. the leave-one-out estimates
3. `theta_J`

The challenge is not the formula. The challenge is deciding what counts as a
valid jackknife deletion for each design/inference family.

## What Can Be Centralized

A generic implementation is plausible at the `InferenceNonParamBootstrap` / `Inference`
layer:

1. compute `theta_hat`
2. compute the leave-one-out distribution
3. return the jackknife-corrected estimate

Conceptually:

```r
compute_jackknife_corrected_estimate = function(unit = c("observation", "cluster", "block", "pair", ...)){
  ...
}
```

The package could also expose:

- `compute_jackknife_bias_estimate()`
- `approximate_jackknife_distribution_beta_hat_T()`
- `compute_jackknife_standard_error()`

The main subclass/design hook would be:

- what the deletion unit is
- how to generate the leave-one-out subset while preserving design semantics

## Main Architectural Gap

For iid fixed designs, delete-1 is obvious.

For this repo, many designs are not iid simple-row settings:

- blocked designs
- paired / matched designs
- KK designs with reservoir structure
- clustered designs
- sequential designs where the observed assignment path matters

So a package-wide jackknife feature needs a policy for deletion units:

- delete 1 subject
- delete 1 block
- delete 1 cluster
- delete 1 matched set
- delete 1 pair / matched pair

This is the main reason jackknife is not “uniformly easy” even though the core
formula is simple.

## Effect On Testing And Confidence Intervals

This is the key inferential distinction.

### Wald

Jackknife correction fits most naturally with Wald inference.

The simplest first implementation would be:

1. replace `theta_hat` by `theta_J`
2. pair it with either:
   - the ordinary standard error, or
   - a jackknife standard error derived from the leave-one-out values

The standard jackknife variance estimate is

```text
Var_hat_J(theta_hat) = (n - 1) / n * sum_i (theta_hat_{(-i)} - theta_bar_(.))^2.
```

So corrected Wald inference is straightforward in principle:

```text
Z_J = (theta_J - delta) / se_J
```

This is the cleanest inferential companion to jackknife correction.

### Score

Score tests are defined by the score under the null model, not by the point
estimate alone.

A jackknife-corrected estimate does **not** define a new score equation.

So unless a separate jackknife-based score-test theory is implemented:

- ordinary score p-values remain unchanged
- ordinary score-inverted confidence intervals remain unchanged

In other words, jackknife correction does **not** naturally update the score
path.

### Likelihood Ratio

Likelihood-ratio tests are defined by the likelihood and the constrained /
unconstrained fits:

```text
LR(delta) = 2 { l(theta_hat) - l(theta_hat_delta) }.
```

Jackknife correction does not change the likelihood or constrained fitting
problem. So it does not automatically define a new LR statistic.

Thus:

- LR p-values would ordinarily remain unchanged
- LR CIs would ordinarily remain unchanged

unless the package chose to build a separate delete-1 calibrated LR procedure,
which would be a different method.

### Gradient

The same logic applies to the gradient test in the package:

- it is defined from the ordinary restricted score and ordinary likelihood-based
  estimate
- jackknife-correcting the point estimate does not automatically define a new
  gradient test

So gradient p-values and CIs would ordinarily remain unchanged.

### Practical Consequence

Like Cordeiro-McCullagh bias correction, jackknife can create a mismatch:

- corrected point estimate
- corrected Wald inference
- uncorrected score / LR / gradient inference

That mismatch is acceptable if the API is explicit, but it should not be hidden.

## Recommended Inference Policy

The cleanest first implementation would be:

1. expose `compute_jackknife_corrected_estimate()`
2. expose jackknife-based Wald inference:
   - `compute_jackknife_wald_two_sided_pval()`
   - `compute_jackknife_wald_confidence_interval()`
3. leave score / gradient / LR methods unchanged unless a path-specific theory
   is added

This is much cleaner than pretending jackknife automatically propagates into all
four asymptotic testing paths.

A practical implementation consequence is that the jackknife standard error
should be computed at the same time as the jackknife-corrected estimate and then
cached on the object. In other words, the package will need a dedicated
jackknife Wald path, not just a point-estimate wrapper.

So the natural implementation pattern is:

1. compute the full-sample estimate
2. compute the leave-one-out estimates
3. form the jackknife-corrected estimate
4. form the jackknife variance / standard error
5. cache both the corrected estimate and its jackknife SE
6. let the dedicated jackknife Wald p-value / CI methods read from that cache

That is preferable to recomputing the leave-one-out distribution separately in
the p-value and CI methods, and it also keeps the estimate / SE pairing
internally coherent.

## Audit Table

The table below is grouped by **computational pattern / design structure**
rather than every alias class separately, because that is what determines
jackknife feasibility here.

### Simple Fixed-Sample Estimators

| Inference path group | Concrete examples | Audit result | Why |
|---|---|---|---|
| Simple mean/risk differences | `InferenceAllSimpleMeanDiff`, `InferenceAllSimpleMeanDiffPooledVar`, `InferenceIncidRiskDiff`, `InferenceIncidNewcombeRiskDiff`, `InferenceIncidMiettinenNurminenRiskDiff` | **Easy** | Leave-one-out recomputation is conceptually direct in simple fixed-sample settings. |
| Simple regression-style asymptotics | `InferenceContinOLS`, `InferenceContinLin`, `InferenceContinRobustRegr`, `InferenceContinQuantileRegr`, `InferenceIncidenceWald` | **Easy** | These already behave like ordinary estimators on subsets, though quantile regression is less smooth. |
| Simple ordinal/rank procedures | `InferenceOrdinalRidit`, `InferenceAllSimpleWilcox`, `InferenceOrdinalPairedSignTest`, `InferenceOrdinalJonckheereTerpstraTest` | **Borderline** | Nonsmooth statistics can make delete-1 estimates noisy, but the mechanics are simple. |
| Exact binomial/Fisher-type procedures | `InferenceIncidenceExactBinomial`, `InferenceIncidExactFisher`, `InferenceIncidExactZhang` | **Borderline** | The corrected estimate is mechanically possible, but the main inferential identity of these methods is exact testing, not bias-corrected point estimation. |

### G-Computation / Marginalization Paths

| Inference path group | Concrete examples | Audit result | Why |
|---|---|---|---|
| G-computation incidence/ordinal/proportion | `InferenceIncidGCompRiskDiff`, `InferenceIncidGCompRiskRatio`, `InferenceOrdinalGCompMeanDiff`, `InferencePropGCompMeanDiff`, KK G-comp variants | **Easy** | The delete-1 estimate is well defined and the package already knows how to refit these objects on subsets. Computational cost is higher but manageable. |

### Fixed-Effects Likelihood GLM-Type Paths

| Inference path group | Concrete examples | Audit result | Why |
|---|---|---|---|
| Binary/count/proportion GLM-like likelihood paths | `InferenceIncidLogRegr`, `InferenceIncidProbitRegr`, `InferenceIncidModifiedPoisson`, `InferenceIncidKKModifiedPoisson`, `InferenceCountPoisson`, `InferenceCountRobustPoisson`, `InferenceCountQuasiPoisson`, `InferenceIncidLogBinomial`, `InferenceIncidBinomialIdentityRiskDiff`, `InferencePropFractionalLogit` | **Easy** | Leave-one-out refits are natural and already aligned with the object model. |
| More bespoke fixed-effects likelihood paths | `InferenceCountNegBin`, `InferencePropBetaRegr`, `InferenceSurvivalWeibullRegr`, fixed-effects ordinal models | **Borderline** | The correction is feasible, but refits are more expensive and small leave-one-out subsets can destabilize optimization. |
| Mixture likelihood paths | zero-inflated / hurdle count, zero/one-inflated beta | **Borderline** | The delete-1 estimate is still definable, but instability and nonconvergence will be more common. |

### Likelihood Paths With Score / Gradient / LR Support

| Inference path group | Concrete examples | Audit result | Why |
|---|---|---|---|
| Standard one-likelihood paths | `InferenceIncidLogRegr`, `InferenceIncidProbitRegr`, `InferenceCountPoisson`, `InferenceContinKKOLSOneLik`, etc. | **Easy** for estimate correction, **Borderline** for inference integration | The estimate correction is easy, but score / gradient / LR methods do not update automatically. |
| Complex likelihood paths | GLMMs, frailty, copula, dependent-censoring transform, custom combined likelihoods | **Difficult** for full integration | Repeated leave-one-out refits are expensive and unstable, and non-Wald testing remains uncorrected anyway. |

### Matching / KK / Structured-Dependence Paths

| Inference path group | Concrete examples | Audit result | Why |
|---|---|---|---|
| KK IVWC / GEE / pass-through estimators | `InferenceAllKKMeanDiffIVWC`, `InferenceIncidKKGEE`, `InferenceCountPoissonKKGEE`, `InferencePropKKGEE`, `InferenceKKPassThrough` descendants | **Borderline** | Delete-1 subject is not obviously the right unit; delete-pair, delete-set, or delete-reservoir-unit may be better. |
| KK one-likelihood / combined-likelihood estimators | `InferenceIncidKKClogitOneLik`, `InferenceCountKKCPoissonOneLik`, `InferenceSurvivalKKLWACoxOneLik`, `InferenceSurvivalKKStratCoxOneLik` | **Borderline** | Mechanically feasible, but deletion semantics and runtime are both nontrivial. |
| Matching / pair-structured exact-like paths | CMH, Extended Robins, matched-set estimators | **Borderline** | These likely need a block/pair jackknife rather than naive subject deletion. |

### Random-Effects / GLMM / Frailty / Copula Paths

| Inference path group | Concrete examples | Audit result | Why |
|---|---|---|---|
| GLMM families | `InferenceIncidKKGLMM`, `InferenceCountKKGLMM`, `InferenceContinKKGLMM`, `InferenceOrdinalKKGLMM`, `InferencePropKKGLMM`, `InferenceOrdinalKKCLMM*` | **Difficult** | Repeated leave-one-out mixed-model refits will be expensive and potentially unstable; the correct deletion unit may also be cluster-level. |
| Frailty / copula survival families | `InferenceSurvivalKKWeibullFrailty*`, `InferenceSurvivalKKClaytonCopula*` | **Difficult** | Same issue with even more bespoke optimization geometry. |

### Survival Nonparametric / Semiparametric Paths

| Inference path group | Concrete examples | Audit result | Why |
|---|---|---|---|
| KM / RMST / log-rank / Gehan-Wilcox | `InferenceSurvivalKMDiff`, `InferenceSurvivalRestrictedMeanDiff`, `InferenceSurvivalLogRank`, `InferenceSurvivalGehanWilcox` | **Borderline** | The delete-1 estimate is natural, but censoring and ties make jackknife variability noisy. |
| Cox / stratified Cox likelihood paths | `InferenceSurvivalCoxPHRegr`, `InferenceSurvivalStratCoxPHRegr`, KK Cox one-likelihood paths | **Borderline** | The corrected estimate is feasible, but score / LR / gradient inference do not inherit a natural jackknife correction. |

### Custom Extension Paths

| Inference path group | Concrete examples | Audit result | Why |
|---|---|---|---|
| User-defined extensions | `InferenceCustomAsymp`, `InferenceCustomRand`, `InferenceCustomBoot` | **Borderline** | Feasibility depends entirely on whether the custom estimator is stable under subset deletion. |

## Easy Tier

The best first targets are:

1. simple fixed-sample mean/risk-difference estimators
2. OLS / Lin / robust regression / standard asymptotic regression estimators
3. g-computation estimators
4. standard fixed-effects GLM-like likelihood paths

These are easiest because:

- the deletion unit is usually just one subject
- subset recomputation is already compatible with the package architecture
- corrected Wald inference can be layered on directly

## Borderline Tier

These are realistic but need policy decisions:

- rank / quantile / nonsmooth estimators
- exact procedures
- negative binomial / beta / Weibull / ordinal likelihoods
- Cox / semiparametric survival
- KK / matching / blocked / clustered structured estimators

The recurring issue is not the jackknife formula. It is:

- what the correct deletion unit is
- how stable the estimator is under deletion
- whether the runtime is acceptable

## Difficult Tier

These are poor targets for an early common rollout:

- GLMMs
- frailty models
- copula models
- dependent-censoring transform models
- the most bespoke combined-likelihood hybrids

For these, the repeated leave-one-out refits are both computationally expensive
and numerically fragile.

## Recommended Implementation Plan

### Phase 1: Explicit Estimate + Wald API

Add explicit methods:

1. `compute_jackknife_corrected_estimate()`
2. `compute_jackknife_bias_estimate()`
3. `compute_jackknife_standard_error()`
4. `compute_jackknife_wald_two_sided_pval()`
5. `compute_jackknife_wald_confidence_interval()`

This is the cleanest starting point because it does not pretend that score / LR
/ gradient paths are automatically corrected.

This phase should also include a dedicated cache for the jackknife estimate-side
quantities, so that one call can populate:

- the ordinary estimate
- the jackknife-corrected estimate
- the jackknife bias estimate
- the jackknife standard error

and subsequent jackknife Wald methods can reuse them directly.

### Phase 2: Standard Delete-1 Families

Implement first for:

1. simple fixed-sample estimators
2. g-computation estimators
3. ordinary fixed-effects GLM-like likelihood paths
4. OLS / Lin / robust regression

These are the most likely to work immediately with the existing
`approximate_jackknife_distribution_beta_hat_T()` helper.

### Phase 3: Structured Deletion Policies

Add optional policies for designs where delete-1 subject is not the right unit:

1. delete-block
2. delete-cluster
3. delete-pair
4. delete-matched-set

Without this phase, package-wide “jackknife correction” would be too naive for
many KK and matching designs.

### Phase 4: Leave Non-Wald Testing Alone At First

The cleanest initial policy is:

- corrected estimate
- jackknife-based Wald inference
- unchanged score / LR / gradient inference

That keeps semantics honest.

## Bottom Line

Jackknife correction is one of the most feasible broad corrections for this repo
because it extends beyond likelihood models and because the package already
contains a shared leave-one-out helper in the bootstrap layer.

The most natural interpretation is:

- a corrected point estimate
- a jackknife-based Wald inference path

Operationally, that means the package should add a dedicated jackknife Wald test
function family and cache the jackknife SE during jackknife estimate
computation.

It is **not** a method that automatically propagates into:

- score tests
- likelihood-ratio tests
- gradient tests

So the right implementation strategy is:

- broad estimate-level rollout
- corrected Wald inference where useful
- selective structured-deletion policies for blocked / matched / clustered
  designs
- no pretense that score / LR / gradient have been jackknife-corrected unless
  separate theory is added
