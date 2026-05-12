# Bayesian Bootstrap Report

## Scope

This report covers how to add a **Bayesian bootstrap** across the package's
bootstrap-capable inference paths by introducing a new class
`InferenceNonParamBayesianBootstrap` that inherits
`InferenceNonParamBootstrap`.

The goal is not to replace the existing nonparametric bootstrap. The goal is to
add a second resampling engine that:

- keeps the observed rows fixed
- draws Bayesian-bootstrap weights instead of duplicated row indices
- reuses as much of the existing p-value / confidence-interval machinery as
  possible

The report focuses on implementation fit with the current codebase.

## Short Answer

This is feasible, but **not** as a pure `bootstrap_sample_indices()` override.

The current bootstrap abstraction in
`EDI/R/inference_all_abstract_non_param_boot.R` is built around:

- drawing bootstrap indices
- physically subsetting / duplicating rows
- recomputing the estimate on the duplicated dataset

A Bayesian bootstrap does not naturally produce duplicated rows. It produces
continuous positive weights, usually Dirichlet-distributed, on the original
observational units or on higher-level resampling units such as pairs, blocks,
or clusters.

So the clean design is:

1. add `InferenceNonParamBayesianBootstrap` as a subclass of
   `InferenceNonParamBootstrap`
2. add one new shared "bootstrap draw spec" abstraction to the base class
3. let the classical bootstrap keep using index-based specs
4. let the Bayesian bootstrap use weight-based specs

That preserves the public bootstrap API and lets `SimulationFramework` continue
working without special casing.

## What The Current Bootstrap Layer Assumes

The central class is `InferenceNonParamBootstrap` in
`EDI/R/inference_all_abstract_non_param_boot.R`.

Today its inner loop is fundamentally:

```r
boot_draw = private$bootstrap_sample_indices(private$n, bootstrap_type)
sub_inf = private$bootstrap_subset_inference(boot_draw, smooth = FALSE)
theta = sub_inf$compute_estimate(...)
```

The important implementation points are:

- `bootstrap_sample_indices()` returns `list(i_b, m_vec_b)` or a bare index
  vector
- `bootstrap_subset_inference()` duplicates the design and inference objects and
  then subsets `X`, `y`, `w`, `dead`, and `m`
- reusable bootstrap workers also assume they can load a bootstrap sample by
  index
- jackknife, BCa, percentile, bootstrap-t, and calibrated bootstrap all build
  on that same index-resample path

This works well for:

- ordinary subject-level bootstrap
- design-aware block bootstrap
- matching bootstrap that resamples matched sets as units and then expands to
  row indices

It is not the right representation for a Bayesian bootstrap.

## Why Bayesian Bootstrap Needs A Different Internal Object

For the Bayesian bootstrap, one bootstrap replicate is:

- the original data
- plus a random weight vector `omega`
- where `omega` is typically Dirichlet on the resampling units

That means the replicate should usually be represented as:

```r
list(
  kind = "weights",
  unit = "subject" | "pair" | "block" | "cluster",
  weights = ...,
  row_weights = ...
)
```

not as:

```r
list(i_b = ..., m_vec_b = ...)
```

Trying to force Bayesian bootstrap weights into duplicated indices would be the
wrong splice for three reasons:

1. Dirichlet weights are continuous, not integer multiplicities.
2. Pair/block/cluster structure should remain explicit instead of being
   approximated by repeated rows.
3. Several estimators could evaluate weighted score equations directly, which is
   cleaner and often faster than materializing a large pseudo-sample.

## Recommended Architecture

### 1. Add A Shared Bootstrap-Draw Spec

Add one new private hook in `InferenceNonParamBootstrap`, conceptually:

```r
bootstrap_sample_spec = function(n, bootstrap_type = NULL){
  list(kind = "indices", i_b = ..., m_vec_b = ...)
}
```

and refactor the existing code so that:

- the current classical bootstrap returns `kind = "indices"`
- the new Bayesian bootstrap returns `kind = "weights"`

Then centralize replication evaluation behind a second hook:

```r
bootstrap_replication_stats_from_spec = function(spec, smooth = FALSE, require_se = FALSE)
```

with default behavior:

- if `spec$kind == "indices"`, use the existing `bootstrap_subset_inference()`
- if `spec$kind == "weights"`, use a new weighted-evaluation path

This is the key splice. Without it, the new subclass will fight the existing
 control flow in multiple places.

### 2. Add `InferenceNonParamBayesianBootstrap`

Create:

```r
InferenceNonParamBayesianBootstrap = R6::R6Class(
  "InferenceNonParamBayesianBootstrap",
  inherit = InferenceNonParamBootstrap,
  ...
)
```

This subclass should override:

- `bootstrap_sample_spec()`
- any bootstrap-type validation needed for Bayesian-bootstrap-specific unit
  choices
- optionally `approximate_jackknife_distribution_beta_hat_T()` if BCa remains
  enabled and a weighted jackknife needs custom handling

It should not need to reimplement the outer CI / p-value methods unless some of
 the more advanced interval types are intentionally restricted.

### 3. Add A Weighted Worker Path

The base class currently has reusable-worker helpers that load a bootstrap
sample by row indices. That should be generalized to accept either:

- index-based bootstrap loads
- weight-based bootstrap loads

Concretely, the reusable worker API likely needs to evolve from:

```r
load_bootstrap_sample_into_worker(worker_state, indices)
```

to something like:

```r
load_bootstrap_spec_into_worker(worker_state, spec)
```

with:

- the old method retained as the `kind = "indices"` implementation
- a new weighted branch that installs row/unit weights into the worker

This keeps the parallel bootstrap paths working for both Unix and Windows.

## The Missing Package-Wide Hook: Weighted Estimation

The real implementation question is not the subclass itself. It is whether each
inference path can evaluate its estimator with externally supplied bootstrap
weights.

That hook does **not** exist today as a generic inference API.

The package already has some weighted fitting primitives:

- `fast_logistic_regression_weighted_cpp`
- `get_logistic_regression_weighted_score_cpp`
- `get_logistic_regression_weighted_hessian_cpp`
- `fast_poisson_regression_weighted_cpp`
- `get_poisson_regression_weighted_score_cpp`
- `get_poisson_regression_weighted_hessian_cpp`

So the infrastructure is partially there. But there is no base-class contract
like:

```r
compute_estimate_weighted = function(row_weights, estimate_only = FALSE)
```

That is the second major splice needed for a package-wide rollout.

## Difficulty Of Threading Observation Weights Through All Inference Paths

This deserves its own audit, because in practice it is the main cost driver.

Adding an observations' weight vector package-wide means more than accepting an
extra argument. For each inference path, all of the following may need to be
weight-aware:

- point-estimate fitting
- constrained/null fitting used by tests or interval inversion
- score and Hessian calculations
- cached design matrices and worker-state reload logic
- randomization/bootstrap helper paths that recompute the estimator internally

For a **first Bayesian-bootstrap implementation**, only the weighted
point-estimation path is essential. Weighted model-based or sandwich standard
errors only become necessary if we later insist on studentized / bootstrap-t
style Bayesian-bootstrap procedures or other methods that require replicate
`se*` values.

So the question is not "can we store `row_weights` somewhere?" The question is
"does the full inferential pipeline remain coherent once row weights are
present?"

### Easy

- logistic-regression incidence paths
- Poisson-regression count/incidence paths
- simple estimators that are already explicit weighted means or weighted score
  roots

Why:

- the repo already has weighted logistic and Poisson C++ fitters, scores, and
  Hessians
- these paths tend to have a single clean model-fitting core
- first-wave Bayesian bootstrap only needs weighted re-estimation, which is
  conceptually straightforward here
- if studentized methods are later desired, these are also the best candidates
  for second-wave support

### Moderate

- log-binomial
- negative binomial
- beta and zero/one-inflated beta
- ordinal likelihood families
- Weibull regression
- g-computation paths

Why:

- these are plausible targets, but many do not currently expose a fully
  parallel weighted fitter / score / Hessian stack the way logistic and Poisson
  do
- some will need new weighted optimizer wrappers, not just light R-side wiring
- some weighted point estimates may be easy while replicate-level standard-error
  support would be materially harder

### Hard

- Cox and stratified Cox paths
- frailty / GLMM / copula paths
- hybrid KK one-likelihood paths
- IVWC / combined estimators that pool multiple sub-estimators
- estimators with bespoke cluster-robust covariance code that assumes unweighted
  contributions

Why:

- these are not just "fit one weighted regression"
- some are based on partial likelihoods, random effects, or custom hybrids
- many would still need bespoke weighted point-estimation implementations
- some pooled estimators would need every constituent sub-estimator to become
  weight-aware first

### Very Hard / Not Worth For Phase 1

- exact inference paths
- heavily custom randomization-based estimators whose estimand is tied to
  permutation structure rather than a weighted likelihood or M-estimation form

Why:

- a Bayesian bootstrap on observational weights may not line up naturally with
  the path's inferential target
- the engineering work would be high and the statistical interpretation may be
  less clean

### Practical Bottom Line

Adding observation weights to **all** inference paths is a large cross-cutting
project. The subclass itself is cheap; making every inference path weight-aware
is the expensive part.

If I had to score the difficulty:

- adding `InferenceNonParamBayesianBootstrap` framework support alone:
  **moderate**
- adding observation-weight support to the easiest first-wave inference paths:
  **moderate**
- adding observation-weight support coherently to all inference paths in the
  repo:
  **hard**

So the implementation plan should explicitly separate:

1. bootstrap framework refactor
2. first-wave weighted inference paths
3. broader package-wide weighted inference adoption

## Recommended Inference-Level API

Add a new private hook in the common inference layer:

```r
supports_weighted_bootstrap = function(){
  FALSE
}

compute_estimate_with_bootstrap_weights = function(row_weights, estimate_only = FALSE){
  stop(class(self)[1], " does not implement weighted bootstrap estimation.")
}
```

Then:

- `InferenceNonParamBayesianBootstrap` calls this hook when the bootstrap spec
  is weight-based
- classes that can implement weighted estimating equations opt in
- classes that cannot opt in yet either:
  - fall back to the existing resampling bootstrap, or
  - reject Bayesian bootstrap explicitly

This lets the package support Bayesian bootstrap broadly without pretending all
paths are equally ready on day one.

For non-likelihood inference paths, the same plan should apply even if there is
no weighted score/Hessian interpretation. The weighted hook should be defined in
terms of "recompute the estimator under bootstrap weights," not "fit a weighted
likelihood."

That means the conceptual contract should be:

```r
compute_estimate_with_bootstrap_weights = function(row_weights, estimate_only = FALSE){
  stop(class(self)[1], " does not implement weighted bootstrap estimation.")
}
```

where `row_weights` may be used by:

- a weighted likelihood fit
- a weighted estimating equation
- a weighted empirical functional
- or a weighted survival estimator such as Kaplan-Meier

So the hook name should stay generic. It should not be framed as
`compute_weighted_likelihood_estimate(...)`, because that would fit only part
of the package.

## Resampling Unit Must Follow The Design

For this package, Bayesian bootstrap weights should not always be at the raw-row
level. They should mirror the unit that the ordinary design-aware bootstrap
treats as exchangeable.

Recommended mapping:

- ordinary fixed or sequential non-blocking designs:
  subject-level Dirichlet weights
- `DesignFixedBlocking`, `DesignFixedOptimalBlocks`, `DesignSeqOneByOneSPBR`:
  either:
  - subject-within-block Dirichlet weights for the analogue of
    `"within_blocks"`, or
  - block-level Dirichlet weights for the analogue of `"resample_blocks"`
- `DesignFixedBlockedCluster`:
  cluster-level Dirichlet weights within stratum for the analogue of
  `"within_blocks"`
- matching designs such as `DesignFixedBinaryMatch` and `DesignSeqOneByOneKK14`
  descendants:
  pair-level Dirichlet weights for matched units, with singleton reservoir units
  handled as singleton clusters

This is important. If a matched-pair design is Bayesian-bootstrapped at the
individual-row level, the replicate no longer respects the design's natural
dependence structure.

## Concrete Implementation Strategy

### Phase 1: Shared Framework

Refactor `InferenceNonParamBootstrap` so the following methods exist:

- `bootstrap_sample_spec()`
- `bootstrap_replication_stats_from_spec()`
- `load_bootstrap_spec_into_worker()`
- `bootstrap_subset_or_weight_inference()` or equivalent internal branching

Keep all public methods unchanged:

- `approximate_bootstrap_distribution_beta_hat_T()`
- `compute_bootstrap_two_sided_pval()`
- `compute_bootstrap_confidence_interval()`

That avoids breaking existing callers and `SimulationFramework`.

For phase one, `bootstrap_replication_stats_from_spec()` should be allowed to
return only `theta*` for Bayesian-bootstrap specs. The package should not force
every Bayesian-bootstrap-supported path to also provide replicate theoretical
SEs from day one.

### Phase 2: New Bayesian Subclass

Implement `InferenceNonParamBayesianBootstrap` with:

- Dirichlet draws using `stats::rgamma(...)` normalization
- unit-aware weight generation based on the design
- row-weight expansion from unit weights

The simplest first constructor surface is probably no new public arguments at
all. The class can just expose Bayesian-bootstrap behavior internally, and
later a public `bootstrap_type = "bayesian"` or `bootstrap_engine = "bayesian"`
switch can be added if desired.

### Phase 3: Easy Families First

The easiest first wave is the likelihood-backed GLM-style paths that already
have weighted fitting support or are very close to it:

- incidence logistic
- incidence modified Poisson / count Poisson
- possibly log-binomial if the optimizer remains stable under weights

These are attractive because the package already contains weighted logistic and
Poisson C++ primitives.

### Phase 4: Broader Asymptotic Families

After the GLM paths, add weighted-bootstrap hooks to:

- g-computation paths
- ordinal regression paths
- selected survival paths
- KK combined / IVWC paths

But these should be audited individually, because some of them are not simple
single-fit estimators.

### Phase 5: Non-Likelihood And Empirical-Functional Paths

After the weighted-likelihood families are working, extend the same bootstrap
spec / weighted-estimation abstraction to non-likelihood paths such as:

- Kaplan-Meier based estimators
- rank-based estimators
- simple contrast / mean / risk-difference functionals
- other estimators that are defined by empirical functionals rather than model
  likelihoods

The implementation idea here is different from the GLM paths:

- do not try to force them through weighted likelihood code
- instead add direct weighted re-estimation methods for each estimator family

For Kaplan-Meier-type paths specifically, the weighted estimator would be based
on a weighted product-limit construction where:

- event counts become weighted event totals
- risk sets become weighted at-risk totals
- the survival curve is recomputed from those weighted risk/event quantities

Then any treatment-effect functional derived from the KM curves can be
recomputed from the weighted curves inside
`compute_estimate_with_bootstrap_weights(...)`.

The same pattern should extend to other non-likelihood paths:

- if the estimator is an empirical mean or contrast, replace raw sums with
  weighted sums
- if the estimator is rank-based, define the weighted analogue of the rank
  statistic or weighted estimating equation
- if the estimator is a two-stage functional, make each stage weight-aware

This is why the package-wide hook should be "weighted bootstrap estimation"
rather than "weighted likelihood fitting." The Bayesian-bootstrap framework
should accommodate both.

## What Is Easy vs Hard

### Easy

These paths are the best first targets:

- logistic-regression incidence paths
- Poisson-regression count and incidence paths
- simple univariate estimators whose estimate is already an explicit weighted
  average or weighted score root

Reason:

- the weighted fitting primitives already exist
- the estimand remains the same when the observed rows are reweighted
- percentile, basic, and bootstrap-SD-based Wald summaries can be formed from
  weighted replicate estimates with little outer-layer change
- if studentized methods are later desired, these are also the best second-wave
  candidates

### Moderate

- negative binomial
- ordinal likelihood paths
- beta / zero-one-inflated beta paths
- Weibull regression
- g-computation paths

Reason:

- these are plausible, but need family-specific weighted fitting or weighted
  post-fit logic
- some rely on custom optimization wrappers rather than an existing weighted
  backend
- they may still be reasonable for a first Bayesian-bootstrap rollout if only
  weighted point estimates are required

### Hard

- Cox partial-likelihood paths
- frailty / GLMM / copula / hybrid one-likelihood KK paths
- estimators whose variance calculation depends on bespoke cluster sandwich
  machinery that is not currently weight-aware

Reason:

- even if the point estimate is definable under weights, the weighted
  re-estimation path may still be materially bespoke
- some of these are better left on the existing resampling bootstrap in phase
  one

## Confidence Intervals And P-Values

For a **first Bayesian-bootstrap implementation**, most of the outer logic can
be reused if the weighted bootstrap produces:

- a replicate estimate `theta*`

That is enough for:

- empirical equal-tail intervals from the weighted-replicate distribution
- percentile-style intervals
- basic intervals
- Wald-like summaries that use the empirical SD of the weighted replicates as
  the uncertainty estimate

Those methods do **not** require each inference path to produce a theoretical
standard error for each weighted replicate.

Methods that **do** require more infrastructure are:

- studentized / bootstrap-t, which need replicate `se*`
- BCa, if we want to preserve the current jackknife-based acceleration path
- calibrated / double-bootstrap / prepivoted, which are nested and more costly
- smoothed bootstrap, which is not a natural first target here

So the distinction should be explicit:

### First implementation

- weighted re-estimation only
- empirical bootstrap distribution of `theta*`
- percentile / basic / equal-tail / bootstrap-SD-based Wald summaries

### Second iteration

- studentized / bootstrap-t where weighted replicate SEs are available
- BCa after validating weighted jackknife behavior
- calibrated / double-bootstrap / prepivoted if they still seem worth the
  engineering cost
- any other methods that require replicate theoretical variance objects

So the pragmatic rollout is:

1. percentile
2. basic
3. bootstrap-SD-based Wald summaries
4. studentized / bootstrap-t where weighted SE is available
5. BCa after validation
6. advanced nested methods later, if at all

## Interaction With `SimulationFramework`

`SimulationFramework` likely does not need structural changes.

As long as:

- the inference object still exposes the same public bootstrap methods
- failures propagate the same way
- the new subclass remains serial/parallel-safe under duplication and worker
  reuse

then the framework should treat Bayesian-bootstrap inference as just another
inference path.

The only likely framework-facing changes are documentation and test coverage.

## Recommended File-Level Changes

### Core bootstrap layer

- `EDI/R/inference_all_abstract_non_param_boot.R`

Add:

- bootstrap draw specs
- weighted worker loading
- weighted replication evaluation
- selective CI-type gating if some methods are unsupported initially

### New class

- `EDI/R/inference_all_abstract_non_param_bayes_boot.R`

Add:

- `InferenceNonParamBayesianBootstrap`
- Dirichlet/unit-weight generation helpers

### Inference base hooks

Likely in:

- `EDI/R/inference_all_abstract.R`
- or the most appropriate shared asymptotic / likelihood abstract layer

Add:

- `supports_weighted_bootstrap()`
- `compute_estimate_with_bootstrap_weights()`

For non-likelihood paths, additional family-specific helpers may be needed in
their own abstract layers, for example weighted survival / weighted empirical
functional helpers, rather than pushing everything into the likelihood base
classes.

### Concrete inference classes

First-wave implementations:

- `EDI/R/inference_incidence_logit.R`
- `EDI/R/inference_count_poisson.R`
- `EDI/R/inference_incidence_modified_poisson.R`

Second-wave candidates:

- `EDI/R/inference_count_negbin.R`
- ordinal likelihood paths
- selected g-computation paths

Later-wave non-likelihood candidates:

- Kaplan-Meier / survival-function paths
- rank-based survival paths
- simple nonparametric contrast estimators

## Proposed Public API

There are two reasonable public surfaces.

### Option A: Separate classes

Create Bayesian-bootstrap inference classes parallel to existing bootstrap
classes, for example:

- `InferenceIncidLogRegrBayesBoot`
- `InferenceCountPoissonBayesBoot`

Pros:

- explicit
- no ambiguity about which bootstrap engine is being used

Cons:

- class proliferation

### Option B: Shared class plus bootstrap-engine switch

Keep the existing inference classes, but add an internal or public switch such
as:

```r
bootstrap_engine = c("classical", "bayesian")
```

Pros:

- fewer class names
- easier `SimulationFramework` integration

Cons:

- more conditional logic inside each class

Given the package's current architecture, I would start with **Option B
internally** and delay any user-facing API decision until the first wave is
working.

## Testing Plan

At minimum, add:

1. unit tests that the Bayesian bootstrap draw spec:
   - sums weights to one
   - respects the design's exchangeable unit
   - expands unit weights to row weights correctly
2. inference tests on easy GLM paths showing:
   - finite bootstrap distributions
   - stable percentile CIs
   - stable bootstrap-SD-based Wald summaries
   - agreement between serial and parallel execution
3. `SimulationFramework` regression tests with Bayesian-bootstrap inference
   paths in both serial and multi-core modes
4. failure-path tests showing unsupported classes error clearly rather than
   silently falling back

## Recommended Rollout

The implementation should be staged.

### Stage 1

- refactor `InferenceNonParamBootstrap` to support bootstrap specs
- add `InferenceNonParamBayesianBootstrap`
- enable weighted-re-estimation Bayesian bootstrap for logistic and Poisson
  families
- support empirical equal-tail / percentile / basic / bootstrap-SD-based Wald
  summaries

### Stage 2

- add weighted replicate standard-error support where it is natural
- enable studentized / bootstrap-t on the first-wave families
- add `SimulationFramework` coverage

### Stage 3

- extend to moderate-difficulty families one at a time
- validate BCa / weighted jackknife behavior

### Stage 4

- extend to non-likelihood paths with direct weighted re-estimation methods
- start with Kaplan-Meier-style estimators and other simple empirical
  functionals

### Stage 5

- decide whether advanced nested methods are worth supporting for Bayesian
  bootstrap at all

## Bottom Line

The right splice is:

- **yes** to a new `InferenceNonParamBayesianBootstrap` subclass
- **no** to implementing it as only an alternate index sampler

The minimal robust design is:

1. generalize the base bootstrap layer from "index resample" to "bootstrap draw
   spec"
2. add a package-wide weighted-estimation hook for inference classes
3. implement Bayesian bootstrap first as weighted re-estimation plus empirical
   bootstrap summaries for the GLM-style paths that already have weighted
   fitting support
4. treat studentized / replicate-SE-based Bayesian-bootstrap features as a
   second iteration
5. keep harder paths on the classical resampling bootstrap until they are
   audited individually

That gives the package a coherent Bayesian-bootstrap architecture without
breaking the existing bootstrap API or overpromising support for difficult
likelihood families.
