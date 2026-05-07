# Smarter Optimization Starts Across the Package

This report describes how to implement stronger optimization starting values
across the package, with the first target set matching the requested rules:

- logistic regression: start at OLS for `y ~ X`
- Poisson and negative binomial: start at OLS for `log(y + 1) ~ X`
- Weibull AFT: start at OLS for `log(y) ~ X` after dropping censored rows
- ordinal cumulative-link models: start slopes at OLS

The package already has a clean separation between:

- R inference wrappers that build the final design matrix
- Rcpp fitters that actually optimize

That separation is useful here. The main implementation choice is whether
smarter starts should be computed in R, in C++, or both.

## Recommendation

Implement the default smart starts inside the Rcpp fitters, but make the
user-facing policy live on the R inference object.

That gives three benefits:

1. Every code path gets the same behavior: main fits, likelihood-ratio null
   fits, bootstrap fits, and randomization fits.
2. The starts are computed on the exact matrix that is actually optimized after
   QR dropping or column filtering.
3. R wrappers stay simple for the first rollout.

For the more complex multi-component models, add small R helpers later to
compose starts from the simpler building blocks.

It is also feasible to expose a user-facing flag such as
`smart_default = TRUE` at the R level. That flag would choose between:

- `smart_default = TRUE`: use the new model-specific smart starts
- `smart_default = FALSE`: preserve the current default starts, which are
  usually all-zero coefficient vectors or similarly crude built-in defaults

I recommend supporting this flag on the R inference interface, not as the
primary public switch on the low-level Rcpp fitters.

## Current State

The current defaults are uneven:

- `fast_logistic_regression_internal()` starts all free coefficients at zero.
- `fast_poisson_internal()` starts from zero for the coefficient vector in the
  LBFGS/Newton path and from a crude working mean in the IRLS path.
- `fast_neg_bin_internal()` starts only the intercept at `log(mean(y))` and
  gets `theta` from a count moment estimator.
- `fast_weibull_regression_cpp()` already uses OLS on `log(y)`, but it uses all
  rows, including censored observations.
- `fast_ordinal_*_regression_cpp()` starts thresholds from generic link-based
  grids and starts all slopes at zero.

So the improvement is not conceptually hard. Most of it is replacing weak
defaults with better deterministic starts and plumbing those starts through the
few fitters that still have no explicit start argument.

## Shared Design

Add a small shared C++ initializer layer in `EDI/src/_helper_functions.h`
and/or a new `EDI/src/optimization_starts.h`. I would keep the helpers pure and
matrix-based.

Recommended helpers:

```cpp
Eigen::VectorXd ols_start_beta(const Eigen::MatrixXd& X,
                               const Eigen::VectorXd& y);

Eigen::VectorXd ols_start_beta_on_log1p(const Eigen::MatrixXd& X,
                                        const Eigen::VectorXd& y);

struct WeibullStart {
    Eigen::VectorXd beta;
    double log_sigma;
};
WeibullStart weibull_aft_start(const Eigen::MatrixXd& X,
                               const Eigen::VectorXd& y,
                               const Eigen::VectorXd& dead);

struct OrdinalStart {
    Eigen::VectorXd alpha;
    Eigen::VectorXd beta;
};
OrdinalStart ordinal_start_from_ols(const Eigen::MatrixXd& X,
                                    const Eigen::VectorXd& y,
                                    edi_ordinal::Link link);
```

Behavior rules:

- Use `CompleteOrthogonalDecomposition` for the OLS solve so starts still work
  after aggressive QR column dropping.
- If the transformed response is constant or non-finite, fall back to zeros.
- If a model has fixed parameters, compute the smart start on the full vector
  and then overwrite fixed entries afterward.
- If there are too few usable rows for the special start, fall back to the
  current default.

The low-level fitters should still accept explicit `start_beta` or
`start_params` where useful, but the user-facing choice between smart and
legacy defaults should be stored on the R inference object and translated into
the appropriate internal start vector before calling Rcpp.

## Final Implementation Spec

The report should be read with the following interface contract fixed. The rest
of the implementation plan should follow these rules without revisiting them
unless there is a later explicit design change:

- `smart_default` lives on the R inference initializer
- low-level fitters still accept explicit `start_beta` or `start_params`
- explicit starts override `smart_default`
- low-level fitters compute smart defaults only when no explicit start is
  supplied and the caller has chosen the smart-start policy

These four points are the fixed interface contract for this project.

This gives a clean split of responsibilities:

- R inference layer:
  owns the user-facing choice between smart and legacy defaults
- low-level fitter layer:
  owns the actual construction and application of explicit or default start
  vectors
- CI inversion / randomization CI layer:
  owns warm-start reuse across nearby `delta` values
- bootstrap orchestration layer:
  may optionally own warm-start reuse across successive resamples

There should be no ambiguity in precedence:

1. if a caller passes `start_beta` or `start_params`, use it
2. otherwise consult the inference object's `smart_default`
3. if `smart_default = TRUE`, build the model-specific smart default
4. if `smart_default = FALSE`, use the legacy default

This precedence rule is the contract:

- `smart_default` is the user-facing default-start policy
- explicit starts are the low-level override mechanism
- explicit starts always win over `smart_default`

## Implemented Framework Notes

The current implementation now supports the following framework pieces:

- `smart_default` is stored on the R inference object initializer and is
  threaded into the low-level fitter calls for the core logistic, Poisson,
  negative-binomial, Weibull, and ordinal cumulative-link inference classes.
- low-level fitters accept explicit `start_beta` or `start_params`, and those
  explicit starts override both `smart_default` and any internally generated
  smart start.
- low-level start construction is performed on the final post-QR fit matrix
  `X_fit`, not on the unreduced `X_full`.
- likelihood-ratio and score-test null-fit evaluation now supports the
  conceptual interface `fit_null(delta, start = NULL)`, with a wrapper that
  caches the previous null-fit solution and reuses it for nearby `delta`
  values during CI inversion and repeated p-value evaluation.
- reusable bootstrap and randomization worker states preserve object-level fit
  warm starts across repeated refits, so bootstrap and randomization-CI delta
  searches can reuse recent solutions when the model family supports explicit
  starts.

Benchmark tooling for this work now lives at
`EDI/inst/benchmarks/benchmark_smart_starts.R`.

## Feasibility of an R-Level `smart_default` Flag

This is feasible.

The main design question is where to expose it in the R interface:

- on the inference-class initializer
- on `compute_estimate()`

### Preferred Location: Initializer

Example shape:

```r
InferenceCountNegBin$new(
  des_obj,
  model_formula = NULL,
  smart_default = TRUE
)
```

or for non-R6 front-end helpers, whatever constructor path ultimately builds the
inference object should carry the same argument.

Why initializer-level is better:

- the choice becomes part of the fit specification for the object
- cached fits, bootstrap refits, randomization refits, and likelihood-ratio
  null fits can all use the same initialization policy automatically
- the object already stores other fitting controls such as optimization
  behavior, so this is a natural place for another fitting policy
- it avoids having the same object produce different optimization behavior
  depending on which method was called last

Implementation shape:

- add `smart_default = TRUE` to the relevant inference-class `initialize()`
  methods
- store it in a private field, for example `private$smart_default`
- in each `generate_mod()` or shared fit helper, use that field to choose the
  start vector before calling the Rcpp fitter

I do not recommend making a mutable setter part of the initial plan. Changing
`smart_default` after construction creates cache-coherence questions and is not
needed to realize the main benefit.

## Explicit Start Vectors and Warm Starts

The framework should explicitly allow caller-supplied starting values, even if
the main user-facing policy is controlled by `smart_default`.

This matters for two reasons:

- advanced callers and composite models may want to provide a start vector
  directly
- inverted score tests and likelihood-ratio tests for confidence intervals can
  benefit substantially from warm-starting the null fit
- bootstrap refits may also benefit from warm starts when successive resamples
  produce nearby optima

### Confidence-Interval Inversion Use Case

For inverted score-test or likelihood-ratio-test confidence intervals, the
framework repeatedly evaluates constrained or null fits at nearby values of a
target parameter, for example through a sequence of `delta` values.

In that setting, the following optimization makes sense:

- solve the first null fit at `delta_1` normally
- solve the next null fit at `delta_2` starting from the fitted parameter
  vector found at `delta_1`
- because the null solution changes smoothly as `delta` varies, the optimizer
  at `delta_2` will usually need far fewer Newton-Raphson or quasi-Newton
  iterations

This is likely the single biggest optimization gain available for confidence
interval inversion.

The same idea can also help randomization confidence-interval construction when
the randomization CI is built by inverting a test over many nearby `delta`
values on the same observed dataset.

### Required Framework Support

To allow that optimization, the plan should require that low-level fitters
accept explicit starting values.

That means:

- coefficient-only models should accept `start_beta`
- full-parameter models should accept `start_params`
- R-level `fit_null(delta)` closures should be able to pass those values
  through to the fitter

The key distinction is:

- `smart_default` chooses how a fit starts when no explicit start is supplied
- `start_beta` or `start_params` lets the caller override that and warm-start
  from a previous solution

So explicit starts are not just optional API polish. They are required for the
CI-inversion optimization path.

### Recommended `fit_null()` Shape

Current conceptual shape:

```r
fit_null = function(delta) {
  ...
}
```

Recommended shape:

```r
fit_null = function(delta, start = NULL) {
  ...
}
```

where:

- `start = NULL` means use the default start policy for that null fit
- `start` supplied means pass that previous solution into the underlying fitter

For coefficient-only models, `start` can be a beta vector.
For models with nuisance parameters, it should be the full parameter vector on
the constrained scale expected by the fitter.

### Where the Warm Start Logic Should Live

The warm-start logic should live above the low-level fitter and below the CI
search routine:

- the low-level fitter only needs to accept explicit starts
- the CI inversion code should keep the last fitted null solution and pass it
  forward to the next nearby `delta`

That separation is important. The optimizer should not try to infer the
previous fit on its own. The CI search routine already knows the sequence of
`delta` values being explored, so it is the right layer to own the warm-start
policy.

### Interaction with `smart_default`

Warm starts and `smart_default` are compatible.

Recommended precedence:

1. if an explicit `start_beta` or `start_params` is supplied, use it
2. otherwise, if `smart_default = TRUE`, use the smart default
3. otherwise use the legacy default

This rule keeps the framework flexible without making it ambiguous.

### Randomization Confidence Intervals

Warm starts can also speed up randomization CI construction, but the gain
depends on where the repeated fitting occurs.

Most promising case:

- the randomization CI algorithm fixes the observed dataset
- it evaluates many nearby null values `delta`
- each new constrained fit starts from the previous constrained solution

In that case, the same smooth-path logic as score/LR inversion applies, and the
speedup can be substantial.

Less promising case:

- the fitting problem changes because the treatment assignment or permuted data
  change from one step to the next

In that setting, the previous fitted solution may be less close to the new
optimum, so warm starts may still help but the gain is less predictable.

So the report should distinguish two uses:

- warm starts across nearby `delta` values on the same observed data:
  likely very helpful
- warm starts across different permutations or treatment reassignments:
  possibly helpful, but less reliable

The framework should still allow both, because supporting explicit starts costs
little once the plumbing exists.

### Bootstrap Refits

Bootstrap refits can also benefit from warm-start support.

Most promising case:

- successive bootstrap samples are fit with the same model structure
- the fitted solution from one resample is reasonably close to the next
  resample solution
- the model is expensive enough that saving a few Newton or LBFGS iterations is
  worthwhile

This is especially plausible for:

- negative binomial
- Weibull and other survival likelihoods
- ordinal cumulative-link models
- bootstrap workflows with many repeated refits on the same dataset shape

Less promising case:

- very unstable bootstrap samples that move the optimum substantially
- tiny or trivial models where the fit is already extremely cheap

So bootstrap warm starts are likely helpful on average for harder model
families, but the gain is less deterministic than warm starts across nearby
`delta` values in CI inversion.

The framework should therefore support bootstrap warm starts, but the plan
should treat them as a secondary optimization after CI-inversion warm starts.

## Model-Specific Plan

### Logistic

Target fitters:

- `fast_logistic_regression_cpp`
- `fast_logistic_regression_weighted_cpp`
- `fast_logistic_regression_with_var_cpp`
- the shared `fast_logistic_regression_internal()`

Implementation:

- Add an optional `start_beta` argument.
- Default it to `ols_start_beta(X, y)`.
- At the R inference level, if `smart_default = FALSE`, use the current zero
  start instead.
- For fixed-parameter fits, apply the fixed values after the OLS start is built.

This is the easiest case. The IRLS solver already supports an arbitrary starting
point through `beta_free`; it currently just initializes at zero.

### Poisson

Target fitters:

- `fast_poisson_regression_cpp`
- `fast_poisson_regression_weighted_cpp`
- `fast_poisson_regression_with_var_cpp`
- `fast_quasipoisson_regression_with_var_cpp`
- the shared `fast_poisson_internal()`

Implementation:

- Add an optional `start_beta` argument.
- Default it to `ols_start_beta_on_log1p(X, y)`.
- At the R inference level, if `smart_default = FALSE`, preserve the current
  legacy start behavior.
- Use the same start for IRLS, Newton, and LBFGS branches.

This is still straightforward because the fitter already has one common internal
entry point.

### Negative Binomial

Target fitters:

- `fast_neg_bin_cpp`
- `fast_neg_bin_with_var_cpp`
- `fast_neg_bin_internal()`

Implementation:

- Add an optional `start_params` argument of length `p + 1`.
- Set `beta_start` to OLS on `log(y + 1)`.
- Keep the current moment-based `theta_start` unless we later decide to build a
  smarter dispersion initializer.
- At the R inference level, if `smart_default = FALSE`, preserve the current
  intercept-plus-moment legacy start.

Recommended default:

```text
start_params = c(beta_ols_log1p, log(theta_moment))
```

This matches the requested behavior and preserves the one part of the current
start that is already useful: the dispersion guess.

### Weibull AFT

Target fitters:

- `fast_weibull_regression_cpp`
- any R helper that fabricates Weibull starts for composite models

Implementation:

- Keep the existing `start_params` API.
- Change the default start generator so OLS on `log(y)` is fit only on rows
  with `dead == 1`.
- Estimate `log_sigma` from the residual standard deviation on those uncensored
  rows.
- If there are too few uncensored rows to fit OLS, fall back to the current
  all-row start or to zeros.
- At the R inference level, if `smart_default = FALSE`, preserve the current
  all-row OLS default.

This is the only requested model that already has the correct shape of API.
Only the default start logic needs to change.

### Ordinal Cumulative-Link Models

Target fitters:

- `fast_ordinal_regression_cpp`
- `fast_ordinal_probit_regression_cpp`
- `fast_ordinal_cloglog_regression_cpp`
- `fast_ordinal_cauchit_regression_cpp`
- corresponding `*_with_var_cpp` wrappers

Implementation:

- Add an optional full `start_params` argument.
- Build `beta_start` from OLS of numeric `y` on `X`.
- Reconstruct thresholds separately so they remain ordered.
- At the R inference level, if `smart_default = FALSE`, preserve the current
  threshold-grid plus zero-slope default.

The main subtlety is that “ordinal models should start at OLS” can only mean
the slope block. These models also require threshold parameters.

Recommended threshold construction:

1. Compute empirical cumulative probabilities
   `p_k = P(Y <= k)` for `k = 1, ..., K - 1`.
2. Convert them to link-scale thresholds using the link quantile:
   `g(p_k)`, where `g` is `qlogis`, `qnorm`, cloglog inverse-quantile, or
   cauchit quantile.
3. Shift all thresholds by a location summary of `X %*% beta_start`
   such as `mean(X %*% beta_start)` so they stay compatible with the nonzero
   slope start.
4. Enforce strict ordering with a small minimum gap.

So the default ordinal start would be:

```text
start_params = c(alpha_from_margins_and_beta_ols, beta_ols)
```

This is still deterministic and cheap, but much better than
“generic threshold grid + zero slopes”.

## R Wrapper Changes

The first rollout can keep most wrappers unchanged if the fitter defaults become
smart. Still, a few R-side updates are worth making.

### Add `smart_default` to Initializers

The inference-class initializers should expose:

```r
initialize = function(..., smart_default = TRUE)
```

for the model families where smarter defaults are introduced.

That is the cleanest public interface for users.

That matters for:

- cached main fits
- bootstrap and randomization refits
- likelihood-ratio null fits
- composite model helpers

The object should store the choice and reuse it consistently in all downstream
fitting paths.

### Keep Explicit Start Vectors Internal

The exported Rcpp wrappers in `EDI/R/RcppExports.R` do not need a user-facing
`smart_default` flag.

They may still gain explicit `start_beta` or `start_params` arguments where
needed for internal reuse, testing, composite-model fitting, and null-fit
warm-starting during CI inversion or bootstrap refitting, but that is a
separate concern from the user-facing API.

### Keep Starts on the Final Fit Matrix

The `generate_mod()` methods in the inference classes already run through
`fit_with_hardened_qr_column_dropping()`. Smart starts should be computed on the
post-drop `X_fit`, not the original `X_full`.

If starts live in C++, this comes for free.

### Reusable Worker Paths

Several classes cache `best_Xmm_colnames` and then refit quickly during
randomization inference. Once the fitter defaults improve, those paths also
benefit automatically without separate work.

## Composite and Downstream Models

After the first rollout, the same start logic should be reused in composite
models rather than inventing new heuristics.

### Zero-Inflated and Hurdle Count Models

Use:

- conditional count block: Poisson or NB smart start
- zero/hurdle block: logistic smart start on `I(y == 0)` or `I(y > 0)`

This is cleaner than the current intercept-only heuristics.

### Clayton-Weibull and Weibull-Frailty Survival Models

Reuse the Weibull AFT start for the regression and scale pieces, then append a
small grid over the dependence or frailty parameter. The package already has
R-side multi-start helpers in `other_helpers.R`; those should be updated to use
the new uncensored Weibull base start.

### Dependent-Censoring Transformation Models

These already have R-side multi-start logic. They can keep that structure, but
their event and censoring regression blocks should be seeded from the same
shared transformed-response helpers instead of ad hoc zeros where possible.

## Suggested R-Level API Changes

I would standardize the inference-class interfaces as follows:

```r
InferenceIncidLogRegr$new(..., smart_default = TRUE)
InferenceCountPoisson$new(..., smart_default = TRUE)
InferenceCountNegBin$new(..., smart_default = FALSE)
InferenceSurvivalWeibullRegr$new(..., smart_default = FALSE)
InferenceOrdinalPropOddsRegr$new(..., smart_default = TRUE)
```

For consistency:

- `smart_default` governs only the internally generated default
- `TRUE` means use the new smarter starts
- `FALSE` means use the current legacy starts

If desired, `compute_estimate()` could support an optional override:

```r
obj$compute_estimate(estimate_only = FALSE, smart_default = NULL)
```

with `NULL` meaning “use the object-level policy”. That is feasible, but again,
the initializer should remain the primary home for this setting.

## Testing Plan

This change should be tested for both correctness and behavior.

### Unit Tests

Add tests that verify the fitter accepts explicit starts and still converges to
the same solution as before on stable problems.

Specific checks:

- logistic with and without intercept
- Poisson with zeros present
- NB with overdispersion
- Weibull with some censoring
- each ordinal link with at least four categories
- each of those with `smart_default = TRUE` and `smart_default = FALSE`
- null-fit sequences for inverted score/LR confidence intervals with and
  without warm starts
- bootstrap refit sequences with and without warm starts

### Initialization Tests

Add direct tests for the helper outputs:

- OLS start length matches parameter dimension
- Weibull start drops censored rows
- ordinal threshold starts are strictly increasing
- fixed-parameter fits preserve the fixed entries after start generation

### Performance and Robustness Tests

The goal is better optimization, so include tests or benchmarks for:

- fewer iterations on representative data
- fewer non-convergence fallbacks
- unchanged final estimates within tolerance

It is enough to test “not worse than before” initially, then optionally add
benchmarks later.

## Estimated Performance Gain

There is a plausible performance gain from this plan, but it should be
described as an estimate, not a guarantee.

The gain will come mainly from:

- fewer optimizer iterations
- fewer line-search or step-halving backtracks
- fewer failed fits that trigger fallback logic
- more stable bootstrap and randomization refits on harder datasets
- potentially faster bootstrap refit sequences when successive resamples have
  nearby optima
- much faster sequences of nearby null fits during confidence-interval
  inversion when warm starts are used

The gain will not come from cheaper per-iteration work. OLS-based starts add a
small up-front cost, so the improvement only pays off if they reduce the number
of expensive likelihood iterations afterward.

For score-test and likelihood-ratio-test inversion, there is a second and often
larger source of gain: warm-starting `fit_null(delta)` from the previous
nearby `delta` solution. That can dominate the default-start improvement,
because CI inversion may require many constrained refits and those refits vary
smoothly along the search path.

Randomization CI construction can benefit from the same mechanism when it
searches over nearby `delta` values on a fixed observed dataset. That gain is
likely smaller or less stable when trying to reuse starts across changing
permutations, because the fitted solution can move more when the assignment
changes.

Bootstrap refits can benefit for the same general reason: successive resamples
often produce nearby, though not identical, optima. The gain is usually less
predictable than CI-inversion warm starts, but it can still be worthwhile for
expensive model families.

### Expected Magnitude by Family

#### Logistic

Expected gain: small to moderate.

- when the current zero start is already close, gain may be negligible
- when treatment and covariates induce a large shift in baseline risk, OLS
  starts should cut several IRLS iterations

Rough estimate:

- typical wall-clock improvement: `0%` to `15%`
- harder problems: up to `20%` to `30%`

#### Poisson

Expected gain: moderate.

- the current zero start implies fitted means near `1` for all rows, which can
  be poor when counts vary substantially
- OLS on `log(y + 1)` should often reduce IRLS or LBFGS work materially

Rough estimate:

- typical wall-clock improvement: `10%` to `25%`
- harder problems: `25%` to `40%`

#### Negative Binomial

Expected gain: moderate to large.

- the current fit already has a reasonable dispersion start, but the slope start
  is weak
- NB likelihood optimization is more expensive than Poisson, so better starts
  have more room to pay off

Rough estimate:

- typical wall-clock improvement: `15%` to `35%`
- harder problems: `30%` to `50%`

#### Weibull AFT

Expected gain: small to moderate for ordinary fits, moderate for censored data
where the current all-row OLS start is distorted.

- because Weibull already uses an OLS-based start, the incremental gain is
  smaller than for Poisson or NB
- the main improvement is using only uncensored rows for the regression start

Rough estimate:

- typical wall-clock improvement: `5%` to `15%`
- heavier censoring settings: `15%` to `30%`

#### Ordinal Cumulative-Link Models

Expected gain: moderate, with the largest benefit in convergence robustness.

- current starts use threshold grids plus zero slopes
- OLS slope starts should reduce the amount of Newton travel required,
  especially with stronger treatment and covariate effects

Rough estimate:

- typical wall-clock improvement: `10%` to `30%`
- harder problems: `25%` to `45%`

### Package-Level Estimate

Across the package features targeted in this report, a realistic first-pass
expectation is:

- median runtime improvement on affected fits: around `10%` to `25%`
- on the more difficult optimization-based families: often `20%` to `40%`
- on already easy, well-conditioned problems: near `0%`

The bigger practical gain may be robustness rather than raw speed:

- fewer non-converged bootstrap replicates
- fewer fallback fits
- less sensitivity to treatment imbalance or large covariate effects

That robustness gain can be more valuable than mean runtime reduction, because a
single failed fit in a bootstrap or randomization loop can dominate user-facing
runtime and reliability.

### What Would Change My Estimate

This estimate assumes:

- moderate-dimensional covariate matrices
- nontrivial signal in the linear predictor
- a noticeable fraction of current fits starting far from the optimum

The estimate would be smaller if:

- most current fits already converge in very few iterations
- datasets are tiny, so optimization cost is negligible
- QR reduction removes most covariates and leaves a nearly trivial model

The estimate would be larger if:

- the package is frequently used in bootstrap or permutation loops
- datasets have strong effects, high censoring, or overdispersion
- current zero starts are causing occasional optimizer instability

### Recommended Validation

To validate the estimate empirically, benchmark before and after on:

- simple low-signal datasets
- moderate-signal datasets
- high-censoring Weibull datasets
- overdispersed count datasets
- ordinal datasets with 4 to 6 levels and strong covariate effects

The most useful summary metrics are:

- median iterations to convergence
- median wall-clock time per fit
- non-convergence rate
- fallback rate in bootstrap or randomization workflows

## Cost of Computing the OLS Start

The smart starts are not free. They add an up-front least-squares solve before
the main likelihood optimization begins.

Usually this cost should be negligible relative to the optimization it improves,
but not always.

### Why It Is Usually Negligible

For the main target families in this report:

- logistic, Poisson, NB, Weibull, and ordinal fits all involve iterative
  likelihood optimization
- each optimization iteration typically evaluates the score, Hessian, or working
  weighted crossproducts on the full dataset
- a single OLS solve is often cheaper than several extra likelihood iterations

So if the OLS start removes even a modest amount of optimization work, it will
usually pay for itself.

This is especially true for:

- negative binomial
- ordinal cumulative-link models
- bootstrap or randomization workflows where hard fits are repeated many times

### When the OLS Cost Can Matter

There are cases where the start-up cost is more visible:

- very small datasets where the original optimizer already converges almost
  immediately
- nearly trivial models, for example intercept-plus-treatment only
- fits where QR reduction drops most columns and leaves a very low-dimensional
  problem
- workflows doing a large number of extremely cheap fits where every constant
  factor matters

In those cases, a smarter start can improve iteration counts while still giving
little or no wall-clock improvement, because the original fit was already too
cheap for initialization savings to matter.

### Practical Interpretation

The relevant tradeoff is:

- smart start cost: one extra transformed-response OLS solve
- smart start benefit: fewer expensive optimization iterations and fewer failed
  fits

For moderate and hard optimization problems, the benefit should dominate.
For very easy problems, the benefit may be neutral or occasionally slightly
negative in wall-clock time.

That is one reason the proposed `smart_default` flag is useful. It provides an
escape hatch for users who care about the absolute fastest path on easy problems
or who are running a workflow where the old default is already sufficient.

### Expected Magnitude of the OLS Overhead

I would describe the overhead qualitatively as:

- usually negligible for NB, Weibull, and ordinal models
- usually small for logistic and Poisson
- potentially noticeable only when the baseline fit is already extremely cheap

In practice, the overhead is most likely to be visible as a small constant-time
penalty, not as a major asymptotic change.

### How to Benchmark This Specifically

If this tradeoff needs to be quantified directly, benchmark:

- total wall-clock time with `smart_default = TRUE` and `FALSE`
- the time spent in start computation alone
- optimizer iteration counts after initialization

The key question is not whether OLS has a measurable cost. It does. The key
question is whether that cost is smaller than the optimization work it saves.
For the intended model families, that should usually be true, but the report
should explicitly acknowledge that it will not be true in every setting.

## Rollout Order

I would do this in four stages.

1. Shared helpers plus logistic and Poisson.
2. Negative binomial and Weibull.
3. All ordinal cumulative-link fitters.
4. Composite models that should reuse the new starts.

That order follows implementation risk. Logistic and Poisson are the easiest,
Weibull is a local default fix, NB is simple once `start_params` is added, and
ordinal requires the most care because thresholds must be rebuilt consistently.

## Risks

The main risks are not numerical in the simple GLM families. They are:

- ordinal threshold construction producing poor or unordered starts
- edge cases with too few uncensored Weibull events
- explicit start plumbing becoming inconsistent across coefficient-only and
  variance-returning wrappers

Those are manageable if the helper layer is kept small and tested directly.

## Bottom Line

This improvement is feasible with modest code change and should be implemented
centrally in the fitters.

The highest-value first pass is:

- Poisson via OLS on `log(y + 1) ~ X`
- negative binomial via OLS on `log(y + 1) ~ X` plus the current `theta` moment
  start
- Weibull via OLS on `log(y) ~ X` using uncensored rows only
- ordinal via OLS slope starts plus link-specific threshold reconstruction

That approach improves the exact model families you named and also creates a
shared starting-value framework that the rest of the package can reuse.

Exposing `smart_default` is feasible and fits this plan cleanly. It should be
implemented primarily on the R inference initializer, with the object carrying
that policy through all downstream fitting work. If a per-call override is ever
added, it should be secondary to the object-level setting, not the main public
API.
