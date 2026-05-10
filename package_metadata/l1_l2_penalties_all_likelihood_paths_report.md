# L1 / L2 Penalties Across All Likelihood Paths

## Bottom Line

Adding **L2 penalties** to all likelihood-backed paths is **feasible but moderately hard**. It is not one change. It is a cross-cutting retrofit touching:

- the R6 likelihood-spec layer
- the native C++ fitters
- score / Hessian / likelihood-ratio test plumbing
- null-fit warm starts and fixed-parameter refits

Adding **L1 penalties** to all likelihood-backed paths is **hard** and, under the current architecture, closer to a **new optimization/inference project** than a routine extension.

The package currently assumes smooth likelihoods with usable gradients and Hessians. That assumption is built into both the optimizer layer and the score / gradient / likelihood-ratio inversion logic. Ridge fits that assumption reasonably well. Lasso does not.

## Why This Is Not Just “Add `lambda`”

The shared likelihood-testing base is in [EDI/R/inference_all_abstract_asymp_lik.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_asymp_lik.R), and it assumes every likelihood path can provide:

- a full fit
- a null fit with the treatment fixed
- a score at the null fit
- an information matrix / Hessian
- a negative log-likelihood for LR tests

That is fine for smooth unpenalized likelihoods. It becomes structurally different once penalties enter:

- **L2** changes the objective, score, and Hessian, but remains smooth.
- **L1** changes the objective to a non-differentiable one, so the current score / gradient / Hessian-based machinery is no longer a clean fit.

There is also a statistical design question that has to be answered first:

1. Are penalties only for optimization stabilization?
2. Are penalties part of the reported estimator?
3. Are they applied to all coefficients or only nuisance covariates?
4. Is the treatment coefficient exempt from penalty?

Without a clear answer, implementation details are underspecified.

## Current Scope

The package currently has **37** R likelihood-spec entry points (`get_likelihood_test_spec`) and **25** native C++ files that call the shared `optimize_fixed_likelihood(...)` path.

That means “all likelihood paths” spans several distinct families, not one homogeneous regression stack:

- ordinary GLM-style models under [EDI/R/inference_all_abstract_mle_or_KM_for_GLMs.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_mle_or_KM_for_GLMs.R)
- count-likelihood models under [EDI/R/inference_count_likelihood.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_likelihood.R)
- companion-likelihood count models such as [EDI/R/inference_count_quasipoisson.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_quasipoisson.R) and [EDI/R/inference_count_robust_poisson.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_robust_poisson.R)
- ordinal likelihood models with threshold parameters
- survival models with partial likelihood / parametric likelihood
- KK combined-likelihood / GLMM / frailty / copula models
- zero-augmented and hurdle models with multiple parameter blocks

So the right question is not “can the package do penalties?” but “how many distinct likelihood engines have to understand penalties?”

## What Makes L2 Plausible

The native optimization layer is centralized enough that ridge can be threaded through it. The key optimizer hooks are in [EDI/src/_helper_functions.h](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/_helper_functions.h), especially:

- `optimize_likelihood_lbfgs`
- `optimize_likelihood_newton`
- `optimize_fixed_likelihood`
- `FixedParameterFunctor`

Most native fitters already expose:

- objective value
- gradient
- Hessian
- fixed-parameter refits

That is exactly what ridge needs. In the simplest form, each model’s negative log-likelihood would become:

`neg_loglik(beta) + 0.5 * lambda2 * sum(beta_penalized^2)`

and its derivatives would be adjusted accordingly.

For models like:

- [EDI/src/fast_logistic_regression.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_logistic_regression.cpp)
- [EDI/src/fast_poisson_regression.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_poisson_regression.cpp)
- [EDI/src/fast_negbin_regression.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_negbin_regression.cpp)
- [EDI/src/fast_ordinal_regression.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_ordinal_regression.cpp)
- [EDI/src/fast_weibull_regression.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_weibull_regression.cpp)
- [EDI/src/fast_coxph_regression.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_coxph_regression.cpp)

that is a real but straightforward implementation pattern.

## What Makes L1 Hard

The current package is built around smooth objectives. Lasso breaks that assumption in three places:

### 1. The optimizer layer

Current solvers are Newton / LBFGS oriented. They do not provide:

- coordinate descent
- proximal gradient
- orthant-wise quasi-Newton
- active-set logic
- KKT-based stopping rules for nondifferentiable penalties

So L1 is not just “add an absolute value term.” It requires either:

- a new generic optimization backend for nonsmooth objectives, or
- model-specific coordinate-descent implementations

### 2. The testing layer

The likelihood-test machinery in [EDI/R/inference_all_abstract_asymp_lik.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_asymp_lik.R) computes:

- score tests
- gradient tests
- likelihood-ratio tests
- inverted confidence intervals

Those are coded as if the objective is differentiable and the null fit has an ordinary score / Hessian. For lasso:

- the score at zero can be a subgradient set, not a vector
- Hessian-based local inversion is not generally valid at active-set changes
- CI inversion by repeatedly refitting penalized lasso models is much less stable

### 3. The parameter layouts

Several families are not simple `beta` vectors:

- ordinal models include thresholds plus regression coefficients
- negative binomial includes dispersion
- zero-inflated / hurdle models have count and zero-process coefficients
- GLMM / frailty / copula models include variance or dependence parameters

For lasso, you must define exactly which blocks are penalized and which are not. That is a per-family design decision, not a generic implementation detail.

## Hard Parts Even for L2

Ridge is feasible, but still not trivial, because the package has several special cases.

### 1. Treatment should probably be unpenalized

Almost every inference class treats the treatment coefficient as the main estimand. Penalizing it would change both:

- the estimand
- the interpretation of the null-fit inversion machinery

In practice, the likely desired behavior is:

- do **not** penalize treatment
- optionally penalize nuisance covariates
- maybe do not penalize intercepts, thresholds, dispersion, or variance parameters either

That means every fitter needs a **penalty mask**, not just a scalar `lambda`.

### 2. Score / LR tests need a policy

Once ridge is added, there are two coherent but different interpretations:

1. Tests are based on the **penalized** objective.
2. Tests remain based on the **unpenalized** likelihood, while the penalty is only used to stabilize the fitted nuisance parameters.

The current codebase does not separate those notions cleanly. For example, `score`, `information`, and `neg_loglik` are all sourced from the same model-specific spec. If ridge is “real,” those methods must all become penalty-aware. If ridge is only for fitting, then the current likelihood-test semantics need a separate unpenalized evaluation path.

### 3. Null-fit warm starts become more important

Many classes already warm-start null fits. Examples include:

- [EDI/R/inference_count_poisson.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_poisson.R)
- [EDI/R/inference_count_negbin.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_negbin.R)
- [EDI/R/inference_count_zero_augmented_poisson_abstract.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_zero_augmented_poisson_abstract.R)

Ridge will usually help convergence, but CI inversion still requires many constrained refits. Those refits must carry the same penalty configuration and mask logic everywhere.

### 4. Companion-likelihood classes are conceptually different

[EDI/R/inference_count_quasipoisson.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_quasipoisson.R) and [EDI/R/inference_count_robust_poisson.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_robust_poisson.R) report quasi / sandwich inference but use a Poisson companion likelihood for score and LR paths.

If you add penalties here, you must decide whether the penalty belongs to:

- the reported quasi estimator
- the companion Poisson likelihood
- both

That is a design choice, not just code plumbing.

## Where The Work Concentrates

The implementation burden is not uniform.

### Low-to-moderate difficulty

These are the best first candidates for ridge:

- logistic regression
- Poisson regression
- log-binomial
- negative binomial
- ordinary ordinal regression
- Weibull regression
- Cox regression

They already have clean native score / Hessian / null-fit paths.

### Moderate-to-high difficulty

- zero-inflated Poisson / ZINB / hurdle models
- GLMM-based paths
- CLMM / ordinal mixed models
- frailty and copula survival models
- KK combined-likelihood models

These have multiple parameter blocks, latent-variable structure, quadrature, or custom combined objectives.

### Highest difficulty

- all L1 versions of the above
- any attempt to preserve the current score / gradient / LR CI inversion semantics under lasso

This is where the current architecture fights the feature.

## Realistic Effort Estimate

## If the goal is ridge only, with treatment unpenalized

This is a **medium-sized refactor**.

Expected work:

1. Add penalty configuration to the public inference objects.
2. Thread penalty metadata through the R likelihood specs.
3. Extend each native fitter signature with penalty strength and penalty mask.
4. Update each model’s objective, score, and Hessian.
5. Decide whether tests use penalized or unpenalized likelihood quantities.
6. Update tests for null-fit inversion, warm starts, and confidence intervals.

That is realistic.

## If the goal is lasso + ridge everywhere

This is a **large project**.

Expected extra work:

1. Add a nonsmooth optimization backend.
2. Redesign score / gradient test support for nonsmooth penalties, or drop those paths for lasso.
3. Define penalty masks per family and parameter block.
4. Rework CI inversion expectations for active-set changes.
5. Add extensive numerical testing because lasso failures will be family-specific.

That is not a quick package-wide enhancement.

## Recommendation

I would not attempt “L1 and L2 for all likelihood paths” in one step.

I would do this in phases:

### Phase 1

Add **ridge only**, and only for the clean native single-likelihood models:

- logistic
- Poisson
- log-binomial
- negative binomial
- ordinary ordinal
- Weibull
- Cox

Also make treatment unpenalized by default.

### Phase 2

Decide the inference semantics:

- penalized objective for tests, or
- unpenalized tests with penalized nuisance fitting

Do not implement further families before this is explicit.

### Phase 3

Extend ridge to the structured models:

- zero-augmented count models
- GLMMs
- frailty / copula models
- KK combined-likelihood classes

### Phase 4

Only then decide whether lasso is worth supporting. My default recommendation is:

- support **ridge package-wide**
- support **lasso only in a limited subset of simple models**, if at all

## Final Assessment

- **L2 across all likelihood paths:** feasible, but not cheap; expect a moderate refactor with several model-family exceptions and a required design decision about inference semantics.
- **L1 across all likelihood paths:** high difficulty; current optimizer and test-inversion architecture are not a natural fit.
- **Best practical path:** ridge first, treatment unpenalized, simple native models first, and treat lasso as a separate project rather than an extension of the same patch.
