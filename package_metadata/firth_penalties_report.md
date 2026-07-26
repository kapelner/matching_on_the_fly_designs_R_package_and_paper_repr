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

So "all likelihood paths" means many different parameter geometries and not one uniform likelihood engine.

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

Without that decision, "add Firth" is underspecified.

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

The closed-form audit below (see [Closed-Form Firth-Gradient Audit](#closed-form-firth-gradient-audit)) confirms this from the opposite direction: a meaningful subset of paths *do* admit a realistic closed-form gradient, so the generic numerical fallback should only be treated as a last resort for the paths that audit rates **No**.

## Why Not Just Use A Derivative-Free Optimizer?

The obstacle above is specifically about getting a *gradient* of the penalty
for use in L-BFGS. A different idea sidesteps gradients entirely: evaluating
`neg_loglik(theta) + 0.5 * log |I(theta)|` is trivial at any `theta`, since
every path already exposes `observed_information` / `fisher_information` as a
function of the parameters. So why not hand that scalar objective to a
derivative-free optimizer (Nelder-Mead, Powell, COBYLA, ...) and skip deriving
adjusted scores or `S3`-style moment tensors altogether?

That genuinely works for **point estimation** on low-dimensional paths, but it
does not turn "Firth everywhere" into a small project, for reasons specific to
how this package uses a fit once it has one:

1. **Refit volume.** CI inversion runs root-finding over repeated *constrained
   null refits* (see [EDI/R/inference_all_abstract_asymp_lik.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_asymp_lik.R)),
   and the test suite / bootstrap paths call these thousands of times.
   Derivative-free methods typically need roughly an order of magnitude more
   function evaluations than L-BFGS/Newton for comparable precision, and each
   evaluation here is not cheap: building `I(theta)` is an `O(p^2)` matrix
   assembly on its own, or a quadrature integral for the GLMM / frailty /
   copula paths. That multiplies the cost most on the paths that are already
   the most expensive.

2. **Precision, and a bug class this project already hit.** Simplex-type
   methods stall at loose tolerances compared to the near-machine-precision
   convergence Newton/L-BFGS give near an optimum. A prior bug fix in this
   project came from exactly this kind of noise: a stochastic p-value
   function feeding `uniroot` had to be pinned via `f.lower`/`f.upper` to stop
   spurious re-evaluation (see the continuous-CI test fixes). Feeding a
   noisier Firth point estimate into that same root-finding machinery would
   reintroduce that failure mode from a different source.

3. **It only solves point estimation, not inference.** Getting a Firth
   estimate derivative-free still leaves score tests, gradient tests, LR
   tests, and Wald SEs needing an actual score/information *at* that
   estimate. Everything under [The Testing Layer Complicates Things](#the-testing-layer-complicates-things)
   is untouched — you would still have to decide what those tests use.

4. **Worst dimensionality exactly where it would be used.** Simplex / pattern
   search methods degrade sharply above roughly 10-20 parameters. The paths
   where closed-form is genuinely impractical — GLMM, frailty, copula,
   ZINB/hurdle — are also the ones with the most parameters (variance
   components, mixture weights, dispersion). So derivative-free is cheapest
   exactly where closed-form already works well (the simple GLMs, where it
   would be strictly worse than the exact gradient already available), and
   most expensive and least reliable exactly where someone would want to use
   it instead.

5. **Behavior near singular `I(theta)`.** `log |I(theta)|` blows up as
   information approaches singularity — precisely the near-separation region
   Firth exists to fix. That is a region of extreme, ill-conditioned
   curvature, and simplex-based methods are known to collapse or stall there
   rather than navigate it robustly, typically with no diagnostic signal that
   it happened.

So a derivative-free optimizer is a reasonable fallback for a handful of
low-dimensional paths if looser tolerances and slower repeated fits are
acceptable, but it relocates the cost from derivation work to runtime and
precision cost — concentrated on exactly the paths that can least afford it.

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

So even if point estimation were implemented, "all likelihood paths" still includes a large amount of inference-surface work.

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

For Cox specifically, you also have the partial-likelihood question: package-wide "Firth" would need a clear stance on partial-likelihood bias reduction versus full-likelihood Jeffreys-type penalization. (A closed-form Firth-Cox derivation is worked out in the [appendix](#appendix-closed-form-s3-risk-set-moment-derivation-for-cox-ph) below.)

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

### If the goal is Firth only for the simplest GLMs

This is a **moderate project**.

A limited rollout for:

- logistic
- Poisson
- maybe log-binomial

looks realistic if you are willing to implement model-specific formulas and testing behavior.

### If the goal is Firth for all likelihood paths

This is a **large project**.

Expected work would include:

1. Define the package-wide statistical meaning of "Firth" for each family.
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

Support Firth only in the most standard native GLM paths, matching the "Yes" rows of the closed-form audit below:

1. `InferenceIncidLogRegr`
2. `InferenceCountPoisson`
3. `InferenceIncidModifiedPoisson` and the Poisson companion-likelihood paths
4. `InferenceIncidLogBinomial`
5. `InferenceIncidBinomialIdentityRiskDiff`
6. `InferenceContinKKOLSOneLik`

### Phase 2

Define the inference semantics explicitly:

- adjusted score behavior
- LR behavior
- CI inversion policy
- Fisher vs observed information

### Phase 3

Only then evaluate whether extension to the "Borderline" tier is worth it — see the [Concrete Conclusions](#concrete-conclusions) below for the likely-worth-it vs. probably-not-early split:

- fixed-effects ordinal cumulative-link models, including `InferenceOrdinalOrderedProbitRegr`
- Cox / stratified Cox
- `InferenceIncidKKClogitOneLik`
- maybe negative binomial

### Phase 4

Treat mixed, zero-augmented, frailty, copula, and KK combined-likelihood paths as separate projects unless there is a very strong reason to unify them under one Firth framework.

## Final Assessment

- **Firth across all likelihood paths:** hard.
- **Why:** the package exposes scores and Hessians broadly, but generic Firth support needs more than that: it needs a coherent information definition and either adjusted-score derivations or derivatives of the information penalty.
- **Best practical path:** restrict Firth to a small subset of simple GLM families first; do not treat package-wide Firth support as a small extension of the existing likelihood base.

---

# Closed-Form Firth-Gradient Audit

## Scope And Decision Rule

This audit covers the package's likelihood-backed inference paths, grouped by the
**actual likelihood engine** they use rather than by every wrapper class that
calls the same score/Hessian code.

Operationally, I audited the paths that participate in the package's
likelihood-backed testing/inversion framework, i.e. the classes/files that
implement `get_likelihood_test_spec()`. That is the relevant surface for the
question "can LBFGS continue uninterrupted under a Firth penalty gradient?".

This means the audit is complete for the current **likelihood-test paths**. It
does **not** claim to cover every asymptotic model in the package, because some
classes use Wald/GEE/pass-through logic without a likelihood-test spec.

The expanded table below lists **every concrete class** on that
likelihood-backed surface. It excludes:

- abstract base classes
- alias names that point to the same concrete class object
- non-likelihood IVWC/GEE/Wald/pass-through classes without `get_likelihood_test_spec()`

The question audited is:

> Does this likelihood path admit a **closed-form Firth / Jeffreys penalty
> gradient** that is realistic enough that the current **L-BFGS** workflow could
> continue without switching to a derivative-free optimizer?

I use three labels:

- **Yes**: a closed-form gradient is standard or straightforward enough that an
  analytic adjusted score / Jeffreys-penalty gradient is realistic.
- **Borderline**: a closed-form gradient exists in principle, but would require
  a bespoke derivation that is algebraically heavy enough that I would not count
  it as a clean "LBFGS continues uninterrupted" path without dedicated work.
- **No**: the path is mixture-, latent-, quadrature-, copula-, or custom
  combined-likelihood enough that a practical closed-form Firth gradient is not
  a credible generic implementation target.

This is a **practical engineering audit**, not a statement about abstract
mathematical existence for arbitrary symbolic differentiation.

## Audit Table

### Incidence

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceIncidLogRegr` | `fast_logistic_regression_cpp` | **Yes** | Standard Bernoulli-logit GLM; classical Firth setting with realistic analytic adjusted score. |
| `InferenceIncidProbitRegr` | `fast_ordinal_probit_regression_cpp` in the 2-category case | **Yes** | Smooth fixed-effects binary probit likelihood; closed-form Jeffreys/Firth gradient is realistic in the same sense as other simple binary GLMs. |
| `InferenceIncidModifiedPoisson` | `fast_poisson_regression_cpp` | **Yes** | Same Poisson likelihood engine as ordinary Poisson regression. |
| `InferenceIncidKKModifiedPoisson` | `fast_poisson_regression_cpp` | **Yes** | Same Poisson companion likelihood engine as above. |
| `InferenceIncidLogBinomial` | `fast_log_binomial_regression_cpp` | **Yes** | Smooth binomial likelihood with noncanonical link; analytic Firth gradient is still realistic. |
| `InferenceIncidBinomialIdentityRiskDiff` | `fast_identity_binomial_regression_cpp` | **Yes** | Smooth binomial likelihood with fixed dispersion; less pleasant algebra than logit, but still tractable. |
| `InferenceIncidKKClogitOneLik` | stacked conditional-logistic + reservoir-logistic path via `fast_logistic_regression_with_var_cpp` | **Borderline** | Smooth logistic-based combined likelihood, but not the ordinary unconditional logit Firth case. |
| `InferenceIncidKKGLMM` | `fast_logistic_glmm_cpp` | **No** | GH-quadrature integrated random-effects likelihood with variance parameter. |
| `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm_cpp` | **No** | Hybrid of conditional logit and quadrature GLMM pieces with shared coefficients. |

### Count

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceCountPoisson` | `fast_poisson_regression_cpp` | **Yes** | Canonical Poisson GLM with explicit score/Hessian. |
| `InferenceCountRobustPoisson` | `fast_poisson_regression_cpp` | **Yes** | Reported estimator is robust/sandwich, but the likelihood-test path is still plain Poisson. |
| `InferenceCountQuasiPoisson` | `fast_poisson_regression_cpp` | **Yes** | Reported estimator is quasi, but the likelihood-test path is still plain Poisson. |
| `InferenceCountNegBin` | `fast_neg_bin_cpp`, `fast_neg_bin_with_var_cpp` | **Borderline** | Smooth full likelihood, but dispersion-parameter derivatives make the Firth gradient bespoke and polygamma-heavy. |
| `InferenceCountZeroInflatedPoisson` | `fast_zero_augmented_poisson_cpp` | **No** | Mixture likelihood with count and inflation blocks; practical analytic Firth gradient is not a clean target. |
| `InferenceCountHurdlePoisson` | `fast_zero_augmented_poisson_cpp` | **No** | Truncation/mixture structure makes the Jeffreys penalty gradient highly bespoke. |
| `InferenceCountZeroInflatedNegBin` | `fast_zinb_cpp` | **No** | Mixture plus dispersion parameter makes closed-form Firth support impractical. |
| `InferenceCountHurdleNegBin` | `fast_hurdle_negbin_cpp` | **No** | Same issue as above with additional hurdle/truncation structure. |
| `InferenceCountKKGLMM` | `fast_poisson_glmm_cpp` | **No** | Random-effects quadrature likelihood; not a practical closed-form Firth path. |
| `InferenceCountKKHurdlePoissonOneLik` | `fast_hurdle_poisson_glmm_cpp` | **No** | Truncated count + random effects + quadrature is too structurally complex. |
| `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined_with_var_cpp` | **No** | Hybrid conditional-plus-marginal likelihood; combined information penalty is bespoke. |

### Continuous

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceContinKKOLSOneLik` | `fast_ols_with_var_cpp` | **Yes** | Gaussian likelihood has explicit information and simple matrix derivatives. |
| `InferenceContinKKGLMM` | `fast_gaussian_lmm_cpp` | **Borderline** | Gaussian structure helps, but variance-component derivatives mean this is no longer a simple drop-in Firth path. |

### Ordinal

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceOrdinalPropOddsRegr` | `fast_ordinal_regression_cpp` | **Borderline** | Fixed-effects ordinal likelihood with thresholds; analytic path exists but needs dedicated derivation. |
| `InferenceOrdinalOrderedProbitRegr` | `fast_ordinal_probit_regression_cpp` | **Borderline** | Same threshold issue as above, with link-specific derivation for probit. |
| `InferenceOrdinalCauchitRegr` | `fast_ordinal_cauchit_regression_cpp` | **Borderline** | Smooth ordinal likelihood, but link-specific adjusted-score derivation is needed. |
| `InferenceOrdinalCloglogRegr` | `fast_ordinal_cloglog_regression_cpp` | **Borderline** | Same as above for the cloglog link. |
| `InferenceOrdinalKKGLMM` | `fast_ordinal_glmm_cpp` | **No** | Cutpoints plus quadrature plus variance parameter make this too bespoke. |

### Proportion

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferencePropBetaRegr` | `fast_beta_regression_cpp` | **Borderline** | Smooth likelihood, but mean-plus-precision structure makes the Jeffreys adjustment nontrivial. |
| `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta_cpp` | **No** | Three-component mixture; not a realistic clean analytic Firth target. |
| `InferencePropKKGLMM` | `fast_logistic_glmm_cpp` | **No** | Same engine and same quadrature/variance-parameter issue. |

### Survival

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceSurvivalCoxPHRegr` | `fast_coxph_regression_cpp` | **Borderline** | Cox bias reduction is plausible, but risk-set derivatives make this a bespoke implementation rather than a plug-in GLM case. |
| `InferenceSurvivalStratCoxPHRegr` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same Cox issue, with extra stratification structure. |
| `InferenceSurvivalKKLWACoxOneLik` | `fast_coxph_regression_cpp` | **Borderline** | Uses Cox partial likelihood over combined data; analytic bias reduction is plausible but still bespoke. |
| `InferenceSurvivalKKStratCoxOneLik` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same as above with strata-specific risk sets. |
| `InferenceSurvivalWeibullRegr` | `fast_weibull_regression_cpp` | **Borderline** | Smooth parametric survival likelihood, but materially more bespoke than the GLM cases. |
| `InferenceSurvivalDepCensTransformRegr` | `fast_dep_cens_transform_optim_cpp` | **No** | Bespoke transformation likelihood with coupled event/censoring parameter blocks. |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty_cpp` | **No** | Frailty integration and variance parameter make analytic Firth support impractical. |
| `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_clayton_weibull_aft_optim_cpp` | **No** | Copula dependence parameter plus Weibull margins plus combined design structure. |

## Concrete Conclusions

### Clean "yes" paths

These are the paths where I would expect a closed-form Firth gradient to be a
realistic extension that preserves the L-BFGS workflow:

- Bernoulli logit GLM
- Poisson GLM
- log-binomial GLM
- binomial identity-link GLM
- Gaussian linear model

These are the best first targets if the goal is to keep the present optimizer
stack and avoid numerical differentiation of the Jeffreys penalty.

### Borderline paths

These paths are analytically smooth enough that a closed-form Firth gradient is
not impossible, but I would not treat them as "LBFGS just keeps going" without
substantial model-specific derivation work:

- negative-binomial regression
- beta regression
- fixed-effects ordinal cumulative-link models, including ordered probit
- Cox / stratified Cox / LWA Cox
- Weibull regression
- Gaussian LMM
- conditional logistic matched-pair combined likelihood

In practice I would split these into two subgroups:

1. **likely worth it**:
   logit/Poisson-adjacent models, some ordinal models, maybe Cox
2. **probably not worth it early**:
   beta, negative binomial with dispersion, Gaussian LMM

### Clear "no" paths

These paths are too mixture- or latent-structure-heavy for a practical generic
closed-form Firth gradient:

- zero-inflated / hurdle count models
- zero/one-inflated beta
- all quadrature GLMM paths
- hurdle Poisson GLMM combined likelihood
- cPoisson combined likelihood
- clogit + GLMM hybrid
- frailty models
- copula models
- dependent-censoring transformation model

For these, a Firth implementation would either become:

- a bespoke research project per family, or
- a numerical penalty-gradient approximation, which defeats the goal of letting
  L-BFGS continue cleanly.

## Recommendation

If the package wants Firth support while preserving the current L-BFGS-based
architecture, I would limit the first implementation set to:

1. `InferenceIncidLogRegr`
2. `InferenceCountPoisson`
3. `InferenceIncidModifiedPoisson` and the Poisson companion-likelihood paths
4. `InferenceIncidLogBinomial`
5. `InferenceIncidBinomialIdentityRiskDiff`
6. `InferenceContinKKOLSOneLik`

After that, the next tier worth evaluating would be:

7. fixed-effects ordinal cumulative-link models, including `InferenceOrdinalOrderedProbitRegr`
8. Cox / stratified Cox
9. `InferenceIncidKKClogitOneLik`
10. maybe negative binomial

I would not plan on package-wide Firth support across the quadrature, mixture,
copula, frailty, and custom combined-likelihood engines if the requirement is
"closed-form gradient so L-BFGS continues uninterrupted."

## Appendix: Closed-Form S3 Risk-Set Moment Derivation For Cox PH

This appendix works out the specific "bespoke" algebra referenced in the
`InferenceSurvivalCoxPHRegr` / `InferenceSurvivalStratCoxPHRegr` rows above,
so that a future implementation does not have to re-derive it. It follows the
Breslow tie-handling and risk-set-sum structure already implemented in
[EDI/src/fast_coxph_regression.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_coxph_regression.cpp)
(the `S0`/`S1`/`S2` accumulation loop around lines 154–214), and shows exactly
what a third moment tensor `S3` would need to add.

### Existing risk-set moments

For each distinct event time `t_k` with risk set `R_k` and Breslow event count
`d_k`, the code already accumulates, per covariate vector `x_i` and Cox weight
`w_i = exp(x_i'β)`:

```
S0(β,t_k) = Σ_{i∈R_k} w_i                       (scalar)
S1(β,t_k) = Σ_{i∈R_k} w_i x_i                    (p-vector)
S2(β,t_k) = Σ_{i∈R_k} w_i x_i x_i'               (p×p matrix)
```

which is exactly `workspace.S1` / `workspace.S2` in the code. From these the
code forms the risk-set mean `x̄_k = S1/S0` (`workspace.e_z`) and the risk-set
weighted covariance `V_k = S2/S0 - x̄_k x̄_k'` (the `hess` accumulation at line
~206–214). Summed with Breslow weights `d_k`, these give the ordinary score
and observed information already returned as `fisher_information`:

```
U(β)  = Σ_k [ Σ_{j∈D_k} x_j - d_k x̄_k ]
I(β)  = Σ_k d_k V_k(β)
```

### Why this is a cumulant-generating-function recursion

`S0(β,t_k) = Σ_{i∈R_k} exp(x_i'β)` is literally the cumulant generating
function (in `β`) of the discrete distribution that puts weight
`exp(x_i'β)/S0` on each risk-set member's covariate vector `x_i` (a per-risk-set
softmax). That means successive `β`-derivatives of `log S0(β,t_k)` are exactly
the successive cumulants of that distribution:

```
∂  log S0 / ∂β        =  x̄_k          (1st cumulant = mean)
∂² log S0 / ∂β∂β'     =  V_k          (2nd cumulant = covariance)
∂³ log S0 / ∂β∂β∂β    =  M3_k         (3rd cumulant = third central moment)
```

`U(β)` and `I(β)` are Breslow-weighted sums of the 1st and 2nd of these. The
Firth penalty gradient needs the derivative of `I(β)`, so it needs the natural
next term in the same recursion: the 3rd cumulant, `M3_k`.

### The S3 moment and the resulting M3 tensor

Define the third weighted moment tensor over the risk set, the direct
extension of `S1`/`S2`:

```
S3(β,t_k) = Σ_{i∈R_k} w_i · (x_i ⊗ x_i ⊗ x_i)     (p×p×p tensor)
```

i.e. `S3_{jlm} = Σ_{i∈R_k} w_i x_{i,j} x_{i,l} x_{i,m}`, computed with the same
incremental risk-set scan already used for `S1`/`S2` (same loop, one more
nested index, same `O(n·p³)` pass over the sorted event/censoring times).

The 3rd cumulant (= 3rd central moment, since this is a natural exponential
family) works out to:

```
M3_{jlm}(β,t_k) = S3_{jlm}/S0 − x̄_j V_{lm} − x̄_l V_{jm} − x̄_m V_{jl} − x̄_j x̄_l x̄_m
```

with `V_{jl} = V_k(β)_{jl}` and `x̄ = x̄_k(β)` as already computed. Equivalently,
`M3_k = ∂V_k/∂β` component-wise: differentiating the risk-set covariance with
respect to each `β_m` reproduces this same tensor — a consistency check on the
formula above.

### Assembling the Firth-adjusted score for Cox

Summing with the same Breslow event-count weights `d_k` used for `I(β)`:

```
∂I(β)/∂β_m = Σ_k d_k · M3_k(β)[:,:,m]              (p×p matrix, one per m)
```

and the Firth penalty gradient follows the generic trace identity:

```
∂/∂β_m log|I(β)| = tr( I(β)^{-1} · ∂I(β)/∂β_m )
```

giving the adjusted score used for the penalized fit / profile score test:

```
U*(β) = U(β) + 0.5 · [ tr(I⁻¹ Σ_k d_k M3_k[:,:,1]), ... , tr(I⁻¹ Σ_k d_k M3_k[:,:,p]) ]
```

This matches the closed-form correction in Heinze & Schemper (2001,
*Biometrics*, "A solution to the problem of monotone likelihood in Cox
regression"), the basis of the `coxphf` R package, confirming the algebra
above is the established Firth-Cox result rather than a novel derivation.

### Why this stays "Borderline" and not "Yes"

- It is a genuine closed form — no numerical third-derivative differentiation
  of `I(β)` is required, so L-BFGS can in principle continue uninterrupted.
- But it requires a new `S3` accumulator (a `p×p×p` tensor) alongside the
  existing `S0`/`S1`/`S2` accumulators, at `O(p³)` storage/compute per risk
  set versus `O(p²)` for the existing Hessian pass — a real cost and code
  change, not a drop-in reuse of the GLM exponential-family shortcut that
  makes logit/Poisson/OLS "Yes".
- `InferenceSurvivalStratCoxPHRegr` / `InferenceSurvivalKKStratCoxOneLik`
  need the same derivation applied independently within each stratum's risk
  sets, and `InferenceSurvivalKKLWACoxOneLik` needs it applied over the
  combined-likelihood risk sets — mechanically the same recursion, but each
  is its own bespoke wiring into the corresponding fitter rather than a
  single shared implementation.
- Efron tie-handling (as opposed to the Breslow approximation assumed above,
  which matches the current `fast_coxph_regression.cpp` implementation) would
  add further per-tie correction terms to `S1`/`S2`/`S3` that are not derived
  here.
</content>
