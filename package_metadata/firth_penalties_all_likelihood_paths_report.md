# Firth Penalties Across All Likelihood Paths

## Bottom Line

Adding **Firth penalties to all likelihood-backed paths** would be **hard**.

It is materially harder than adding a simple L2 penalty, but easier in principle than adding L1. The reason is:

- Firth correction keeps the objective smooth.
- But it does **not** just add a scalar ridge term.
- It adds a parameter-dependent penalty based on the information matrix, typically
  `-ell(theta) - 0.5 * log |I(theta)|` in neg-loglik form.

So every supported likelihood path would need to understand not only:

- the objective
- the score
- the Hessian

but also the **information matrix as a function of the parameters**, and usually the derivative of the Firth penalty as well.

Under the current architecture, that is a large cross-cutting change rather than a localized retrofit.

## Why Firth Is Different From Ridge

Ridge fits the current optimization setup fairly well because each model can just change:

- objective
- gradient
- Hessian

by adding a simple quadratic term.

Firth is different because the penalty itself depends on the information matrix:

`0.5 * log |I(theta)|`

That creates two immediate complications:

1. You need a clear definition of `I(theta)`.
   Usually this is the expected Fisher information, but some implementations use observed information or model-specific adjusted-score formulas.

2. To optimize the penalized objective, you need the derivative of that penalty.
   The derivative of `log |I(theta)|` depends on derivatives of `I(theta)` with respect to the parameters.

That means first and second derivatives are no longer enough in the generic optimizer layer.

## What The Current Code Supports Well

The package already has a broad likelihood-testing abstraction in [EDI/R/inference_all_abstract_asymp_lik.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_asymp_lik.R). Each likelihood-backed class supplies some combination of:

- `fit_null`
- `score`
- `observed_information`
- `fisher_information`
- `neg_loglik`

There are also shared native optimizer helpers in [EDI/src/_helper_functions.h](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/_helper_functions.h), including:

- `optimize_likelihood_lbfgs`
- `optimize_likelihood_newton`
- `optimize_fixed_likelihood`
- `numerical_gradient`
- `numerical_hessian`

That is enough for ordinary smooth likelihoods and even some simple smooth penalties.

It is **not** enough for generic Firth support.

## The Main Technical Obstacle

For direct optimization of a Firth-penalized objective, the optimizer would need the gradient of:

`neg_loglik(theta) + 0.5 * log |I(theta)|`

The gradient of the second term is not just a function of `I(theta)`. It depends on how `I(theta)` changes with `theta`. In practice that means:

- third derivatives of the log-likelihood, or
- model-specific adjusted-score formulas, or
- repeated numerical differentiation of the information matrix

The current codebase usually stops at:

- objective
- score
- Hessian / information

That is the key reason this is hard.

## Current Scope

The package currently has **37** R likelihood-spec entry points (`get_likelihood_test_spec`) and **25** native C++ files using the shared `optimize_fixed_likelihood(...)` path.

That spans several distinct model families:

- ordinary GLM-style models
- count-likelihood models
- companion-likelihood count models
- ordinal likelihood models
- survival models
- KK combined-likelihood models
- GLMM / frailty / copula paths
- zero-augmented and hurdle models

So “all likelihood paths” means many different parameter geometries and not one uniform likelihood engine.

## Information-Matrix Consistency Is A Real Problem

Canonical Firth correction is tied to Jeffreys-prior penalization, which is usually expressed through the **expected Fisher information**.

But in this codebase, information support is mixed:

- many paths expose only observed information
- a smaller subset explicitly advertises Fisher support
- some paths treat observed information as the main information object

The default likelihood base in [EDI/R/inference_all_abstract_asymp_lik.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_asymp_lik.R) defaults to observed information unless a class explicitly supports Fisher.

Only a small subset of classes explicitly overrides `supports_fisher_information()`. So before implementing Firth package-wide, you would have to decide:

1. Is Firth defined here using expected Fisher information?
2. Is observed information an acceptable substitute in some families?
3. Is there one package-wide rule, or a per-family rule?

Without that decision, “add Firth” is underspecified.

## Why Generic Numerical Firth Is Not Attractive

There is a numerical gradient helper and a numerical Hessian helper in [EDI/src/_helper_functions.h](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/_helper_functions.h).

In theory, you could try to build a generic Firth objective by:

1. evaluating `I(theta)`
2. computing `log |I(theta)|`
3. numerically differentiating that penalty

But for this package, that would be a weak approach:

- it would be slow for CI inversion and null-fit refits
- it would be fragile near singular information matrices
- it would be especially brittle in high-parameter or multi-block models
- it would interact poorly with repeated constrained fits used by score / gradient / LR inversion

So while a numerical fallback might work for prototypes, it is not a good package-wide answer.

## Existing Penalty Support Does Not Solve This

Some native model files already include ad hoc smooth penalties or barriers. For example:

- [EDI/src/fast_logistic_glmm.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_logistic_glmm.cpp)
- [EDI/src/fast_hurdle_poisson_glmm.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_hurdle_poisson_glmm.cpp)

Those are useful evidence that penalized objectives can be carried through the current solver stack.

But they do **not** amount to a general penalty framework:

- the penalties are local to specific models
- the score and Hessian are manually adjusted
- the reported `neg_loglik` is corrected back to the unpenalized quantity afterward

Firth would need a general policy for:

- fitting
- reported log-likelihood
- null-fit inversion
- score tests
- gradient tests
- likelihood-ratio tests

That policy does not exist yet.

## The Testing Layer Complicates Things

The likelihood machinery in [EDI/R/inference_all_abstract_asymp_lik.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_asymp_lik.R) supports:

- Wald tests
- score tests
- gradient tests
- likelihood-ratio tests
- inverted confidence intervals via repeated null refits

Firth complicates all of them.

### 1. Score tests

If the estimator becomes Firth-penalized, what is the score test supposed to use?

- the ordinary score?
- the adjusted score?
- expected Fisher or observed information?

Those are not interchangeable.

### 2. Likelihood-ratio tests

Should LR tests compare:

- penalized likelihoods, or
- unpenalized likelihoods evaluated at Firth estimates?

Those answer different questions.

### 3. CI inversion

The package repeatedly computes constrained null fits. With Firth, every null fit also becomes information-penalized, which means:

- more expensive optimization
- greater sensitivity to near-singular information
- more places where model-specific behavior matters

So even if point estimation were implemented, “all likelihood paths” still includes a large amount of inference-surface work.

## Family-Specific Difficulty

### Lower difficulty

These are the most plausible first targets:

- logistic regression
- Poisson regression
- perhaps log-binomial

These are the classical settings where Firth-style bias reduction is most familiar, and they already have clean native score / Hessian paths.

Even here, though, a robust implementation is not trivial unless you derive or code the adjusted score explicitly.

### Moderate difficulty

- negative binomial
- ordinary ordinal regression
- Weibull regression
- Cox regression

These still have manageable native likelihood code, but Firth is no longer a routine textbook add-on in the same way as logistic regression.

For Cox specifically, you also have the partial-likelihood question: package-wide “Firth” would need a clear stance on partial-likelihood bias reduction versus full-likelihood Jeffreys-type penalization.

### High difficulty

- zero-inflated and hurdle models
- ZINB
- zero/one-inflated beta
- mixed models
- GLMMs
- frailty models
- copula survival models
- KK combined-likelihood models

These models have multiple parameter blocks, latent structure, quadrature, or custom combined objectives. Defining and differentiating the right information penalty for them is much more involved than for ordinary GLMs.

## The Treatment-Effect Target Also Matters

The package is organized around inference on a treatment coefficient. Firth correction is not just a nuisance stabilization term; it changes the estimator itself.

So you have to decide:

1. Is the package reporting the Firth-adjusted treatment estimate as the new estimand?
2. Is Firth only a separation/bias-reduction fallback for nuisance stabilization?
3. Does Firth apply to all parameters, including treatment, thresholds, dispersion, and variance parameters?

Those decisions are especially important in:

- ordinal threshold models
- negative binomial dispersion models
- GLMM variance-parameter models
- zero-augmented models with multiple submodels

This is another reason the feature does not reduce to a generic base-class patch.

## Realistic Implementation Shapes

There are two plausible implementation strategies.

### Strategy 1: Model-specific adjusted-score implementations

This is the statistically clean route.

For each supported family, implement:

- Firth-adjusted score
- suitable information matrix
- corresponding constrained null fits

Pros:

- numerically credible
- fast enough for repeated refits
- matches the literature better

Cons:

- a lot of separate work
- little code reuse across heterogeneous models

### Strategy 2: Generic information-penalized objective

Build a common Firth-style wrapper around the information matrix and optimize:

`neg_loglik(theta) + 0.5 * log |I(theta)|`

Pros:

- conceptually unified

Cons:

- difficult gradient implementation
- likely reliance on finite differences
- fragile near singularities
- expensive during CI inversion

For this codebase, Strategy 2 looks much less attractive than it sounds.

## Practical Effort Estimate

## If the goal is Firth only for the simplest GLMs

This is a **moderate project**.

A limited rollout for:

- logistic
- Poisson
- maybe log-binomial

looks realistic if you are willing to implement model-specific formulas and testing behavior.

## If the goal is Firth for all likelihood paths

This is a **large project**.

Expected work would include:

1. Define the package-wide statistical meaning of “Firth” for each family.
2. Decide on Fisher vs observed information policy.
3. Add Firth-capable fitters model by model.
4. Update null-fit wrappers for score / gradient / LR inversion.
5. Decide whether tests are penalized, adjusted, or unpenalized-at-adjusted-estimates.
6. Add heavy numerical testing for edge cases and singular-information behavior.

That is well beyond a routine refactor.

## Recommendation

I would **not** try to add Firth penalties to all likelihood paths at once.

I would do this in phases:

### Phase 1

Support Firth only in the most standard native GLM paths:

- logistic
- Poisson
- maybe log-binomial

### Phase 2

Define the inference semantics explicitly:

- adjusted score behavior
- LR behavior
- CI inversion policy
- Fisher vs observed information

### Phase 3

Only then evaluate whether extension to:

- negative binomial
- ordinal regression
- Weibull / Cox

is worth it.

### Phase 4

Treat mixed, zero-augmented, frailty, copula, and KK combined-likelihood paths as separate projects unless there is a very strong reason to unify them under one Firth framework.

## Final Assessment

- **Firth across all likelihood paths:** hard.
- **Why:** the package exposes scores and Hessians broadly, but generic Firth support needs more than that: it needs a coherent information definition and either adjusted-score derivations or derivatives of the information penalty.
- **Best practical path:** restrict Firth to a small subset of simple GLM families first; do not treat package-wide Firth support as a small extension of the existing likelihood base.
