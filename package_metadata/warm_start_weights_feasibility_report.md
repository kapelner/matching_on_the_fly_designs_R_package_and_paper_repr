# Warm-Start Weights Feasibility Report

## Scope

This report evaluates the feasibility of adding a `warm_start_weights`
mechanism for likelihood paths.

The proposed idea is:

- instead of warm-starting the optimizer with a previous parameter vector such
  as `start_beta` or `start_params`
- initialize the optimizer using a previous set of **IRLS working weights**
- where these are the algorithmic weights from the iterative solver, not
  externally supplied observation weights

This is a narrower and more technical question than the Bayesian-bootstrap
report. The question here is not whether weights can enter the likelihood. The
question is whether the package should replace or augment parameter warm starts
with **working-weight warm starts**.

## Short Answer

This is only feasible for a limited subset of inference paths, and even there
it is usually **less attractive than parameter warm starts**.

Why:

- IRLS working weights are not free parameters
- they are deterministic functions of the current linear predictor or fitted
  mean
- in the current code, the solvers are organized around warm-starting
  coefficients / parameter vectors, not warm-starting an IRLS state
- most non-GLM likelihood paths do not have IRLS working weights at all

So a `warm_start_weights` feature is:

- potentially useful for a few GLM-style paths
- not a package-wide substitute for `start_beta` / `start_params`
- best viewed as an optional solver optimization for the IRLS-capable families,
  not as a general abstraction for all likelihood inference

## Important Distinction

There are three different objects here:

1. **observation weights**
   User- or bootstrap-supplied case weights that enter the likelihood.
2. **parameter warm starts**
   Previous coefficient / parameter vectors used as optimizer initialization.
3. **IRLS working weights**
   Algorithmic weights recomputed each iteration, usually from the current
   fitted mean.

The proposal in this report concerns item 3 only.

## What The Current Code Does

The existing likelihood paths overwhelmingly warm-start using parameter vectors.

Examples:

- `set_fit_warm_start(...)` and `get_fit_warm_start_for_length(...)` in
  [inference_all_abstract.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract.R:348)
- logistic paths passing `start_beta`; see
  [inference_incidence_logit.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_incidence_logit.R:139)
- negative-binomial paths passing `start_params`; see
  [inference_count_negbin.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_negbin.R:138)

In the IRLS-capable GLM engines, working weights are computed from the current
iterate:

- logistic: `w_diag = mu * (1 - mu)` possibly multiplied by observation weights
  in [fast_logistic_regression.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_logistic_regression.cpp:118)
- Poisson: `w = mu` or `w = mu * weights` in
  [fast_poisson_regression.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_poisson_regression.cpp:117)
- log-binomial: path-specific `working_weights` are carried explicitly in
  [fast_log_binomial_regression.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_log_binomial_regression.cpp:168)

So today the codebase is parameter-warm-start-centric, not IRLS-state-centric.

## Feasibility By Path

The table below separates:

- paths whose model-fitting core truly has IRLS working weights
- paths that only touch IRLS-capable fitters internally
- paths for which `warm_start_weights` is not a natural abstraction

| Path / family | IRLS working weights in optimizer? | Feasibility of `warm_start_weights` | Why |
|---|---|---|---|
| `InferenceIncidLogRegr` | Yes | Moderate | Logistic IRLS weights exist and are cheap to cache, but `start_beta` already determines the initial `mu` and therefore the initial weights. |
| `InferenceCountPoisson` | Yes | Moderate | Same comment; Poisson IRLS weights are derived from `mu`. |
| `InferenceIncidModifiedPoisson` | Yes | Moderate | Same Poisson fitter underneath, so technically feasible but probably incremental relative to `start_beta`. |
| `InferenceCountQuasiPoisson` | Yes for point-estimate fit | Moderate | Point estimate uses Poisson-style GLM fitting, but variance path is quasi-specific. |
| `InferenceCountRobustPoisson` | Yes for point-estimate fit | Moderate | Feasible for the fit, though robust variance still depends on post-fit quantities. |
| `InferencePropFractionalLogit` | Yes | Moderate | Uses logistic fitting, so working-weight initialization is available in principle. |
| `InferenceIncidLogBinomial` | Yes, custom path | Moderate to hard | Has working weights, but this fitter is less standard and often more numerically delicate. |
| `InferenceIncidGComp*` | Indirectly | Moderate | Uses logistic sub-fits, but the overall estimator is more than one IRLS solve. |
| `InferencePropGComp*` | Indirectly | Moderate | Same issue as incidence g-comp. |
| `InferenceIncidKKGComp*` | Indirectly | Moderate | Same issue plus KK structure. |
| `InferenceIncidKKMarginal*` | Indirectly via Poisson fits | Moderate | Feasible for the GLM core only. |
| `InferenceCountKKCPoisson*` | Partly | Hard | Some internal components are IRLS-capable, but the full path is bespoke and not a clean single-IRLS object. |
| Probit / ordinal likelihood paths | Usually no | Low | These are not organized around GLM IRLS working weights in the same way. |
| Negative binomial | No | Low | Uses parameter optimization, not a shared IRLS-state abstraction. |
| Cox / stratified Cox | No | Low | Partial-likelihood optimization, not IRLS working-weight iteration. |
| Weibull / frailty / copula / GLMM paths | No | Low | Custom optimization or quadrature-heavy paths, not IRLS-weight driven. |
| Combined / IVWC KK likelihood paths | Usually no | Low | Multiple sub-fits or pooled estimators make a single `warm_start_weights` abstraction weak. |

## Why `warm_start_weights` Is Usually Redundant With `start_beta`

For the standard GLM-style IRLS families:

- the working weights are determined by the current fitted mean
- the fitted mean is determined by the current parameter vector

So if you already have a good `start_beta`, you automatically get:

- good initial `eta`
- good initial `mu`
- good initial working weights

In other words, for logistic and Poisson regression:

```text
start_beta  ->  eta^(0)  ->  mu^(0)  ->  working_weights^(0)
```

That makes `warm_start_weights` somewhat redundant unless one of the following
is true:

1. you want to warm-start from an approximate previous fit without trusting the
   previous coefficient vector directly
2. the solver is being refactored to carry a richer IRLS state than just
   coefficients
3. there is a path where good working weights matter but coefficient warm starts
   are unavailable or unstable

Those are possible, but they are narrower than "replace warm starts with
working weights."

## What Would Need To Change

To support `warm_start_weights` cleanly, the IRLS-capable fitters would need a
new API surface.

Conceptually:

```r
fast_logistic_regression_cpp(
  X, y,
  start_beta = NULL,
  warm_start_weights = NULL,
  ...
)
```

or internally:

```r
list(
  beta0 = ...,
  mu0 = ...,
  working_weights0 = ...
)
```

This is a meaningful solver change because the current fitters are written to
take:

- a parameter vector
- then derive the working weights internally

They are not currently written to accept an external IRLS state.

## Engineering Risks

### 1. Weight-only warm starts may be underidentified

For logistic or Poisson regression, working weights alone do not determine the
current linear predictor. Two different coefficient vectors can produce similar
working weights.

So if `warm_start_weights` is provided without a matching `start_beta`, the
solver still needs some coherent way to recover or guess:

- initial `eta`
- initial `mu`
- or initial pseudo-response

That makes a pure weight-only warm start weak.

### 2. It complicates the fitter API

Today the fitters mostly accept:

- `start_beta`
- `start_params`
- `smart_start`

Adding `warm_start_weights` means:

- more branching inside the solver
- more validation logic
- decisions about precedence when both `start_beta` and `warm_start_weights`
  are given

### 3. Limited reuse outside GLM-style paths

A package-wide abstraction is only worthwhile if many classes can use it. Here,
most of the likelihood paths would not benefit.

## When It Could Be Worthwhile

There are two plausible use cases.

### 1. Repeated nearby weighted GLM fits

This is the strongest case.

Examples:

- Bayesian-bootstrap GLM replicates
- classical weighted bootstrap experiments
- repeated null-constrained refits in a GLM family

In those settings, carrying forward a richer IRLS state might shave iterations.

### 2. Solver diagnostics / reproducibility

If the package ever wants to expose the internal fit state for debugging or
research, caching:

- `beta`
- `mu`
- `working_weights`

could be useful.

But that is still an argument for **augmenting** parameter warm starts, not
replacing them.

## Recommended Design If Implemented

If this feature is attempted, it should be narrow.

### Recommendation 1

Do **not** replace `start_beta` / `start_params` with `warm_start_weights`.

### Recommendation 2

If implemented, make it an **optional augmentation** for the IRLS-capable GLM
fitters only:

- logistic
- Poisson
- modified Poisson
- quasi/robust Poisson
- fractional logit
- maybe log-binomial after stability testing

### Recommendation 3

Prefer a richer cached fit state instead of raw weight-only starts, for example:

```r
list(
  beta = ...,
  mu = ...,
  working_weights = ...
)
```

Then the solver can:

- use `beta` as the canonical warm start
- optionally skip recomputing some derived quantities when safe

That is more defensible than a bare `warm_start_weights` vector.

## Bottom Line

`warm_start_weights` is not a strong package-wide abstraction.

It is:

- technically feasible for a limited set of GLM-like IRLS paths
- mostly redundant with good parameter warm starts
- not naturally portable to the majority of likelihood families in the repo

So the practical recommendation is:

1. keep `start_beta` / `start_params` as the primary warm-start mechanism
2. if needed, add a richer cached IRLS-state warm start for the small set of
   IRLS-capable GLM families
3. do not design the broader likelihood framework around `warm_start_weights`
   alone
