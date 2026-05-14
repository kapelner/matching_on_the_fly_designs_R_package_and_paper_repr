# R6 MixIn Classes Report

## Short answer

Using MixIn-style composition in this package is plausible for some of the current **behavior-bundle** inference classes, but not for the whole `Inference*` hierarchy.

The current codebase uses inheritance for three different jobs at once:

1. semantic type hierarchy
2. code reuse
3. runtime capability detection via `inherits(...)`, `is(...)`, and `get_inherit()`

MixIns help with job 2, but they weaken jobs 1 and 3 unless we add a new explicit capability/role system.

So the right conclusion is:

- keep the core `Inference -> ...` semantic spine as real inheritance
- consider MixIns only for reusable behavior bundles
- do **not** expect copied-in MixIns to satisfy current `is(...)` / `inherits(...)` checks automatically

## What “MixIn” would mean in this repo

Base `R6` does not give us true multiple inheritance. In practice, “MixIn” here would mean one of these patterns:

1. a helper that returns `public` / `private` method lists to splice into `R6::R6Class(...)`
2. a helper that mutates a generator after definition
3. composition through delegated helper objects

All three are forms of **behavior composition**, not extra superclass nodes in the real R6 inheritance chain.

That distinction matters because the package currently relies heavily on the real inheritance chain.

## Current code that depends on true inheritance

The package does not only use inheritance for code sharing. It also uses it for runtime logic.

Important examples from the current source:

- [EDI/R/inference_suite.R](../EDI/R/inference_suite.R) walks `get_inherit()` to decide whether a generator is an `Inference` subclass and which ones are abstract infrastructure classes.
- [EDI/R/inference_all_abstract_rand_ci.R](../EDI/R/inference_all_abstract_rand_ci.R) uses `inherits(self, "...")` to decide whether an object should be treated as “GLM-like” for transformed randomization CI logic.
- [EDI/R/inference_all_abstract_rand.R](../EDI/R/inference_all_abstract_rand.R) uses `inherits(self, "InferenceAbstractKKQuantileRegrIVWC")` and `inherits(self, "InferenceAbstractKKQuantileRegrOneLik")` for branch-specific behavior.
- Many initializers use `class(self)[1]` for errors and policy dispatch; that part would still work, but only for the concrete class, not for MixIn membership.

Therefore MixIns are not a drop-in substitute for the current inheritance tree.

## Classes that should stay in the true inheritance spine

These classes are serving as actual semantic bases and API-defining layers, not just code-sharing containers:

- `Inference`
- `InferenceRand`
- `InferenceRandCI`
- `InferenceNonParamBootstrap`
- `InferenceAsymp`
- `InferenceAsympLik`
- `InferenceExact`

Why these should stay real superclasses:

- they define broad public API families
- they are part of discovery in `InferenceSuite`
- they are meaningful categories to test with `inherits(...)`
- they describe what an object fundamentally **is**, not just what helper logic it happens to use

`InferenceParamBootstrap` is currently very thin, but it is also a plausible semantic node if parametric LR bootstrap is added later. It can remain a real superclass if we want a discoverable “this family supports or may support parametric null bootstrap” branch.

## Best MixIn candidates

These are the current classes most naturally interpreted as reusable behavior bundles rather than true ontological categories.

   ### 1. `InferenceAsympLikStdModCache`

   This is the strongest MixIn candidate.

   Today it mainly contributes:

   - the standardized `generate_mod()` contract
   - `shared()` caching of `beta_hat_T`, `s_beta_hat_T`, `df`, and `cached_mod`
   - bootstrap worker plumbing for design-backed refits
   - default score / gradient / LR dispatch through `get_likelihood_test_spec()`

   That is a reusable implementation bundle, not a natural semantic type.

   If converted to a MixIn, it would likely become something like:

   - `StdModelCacheMixin`
   - `DesignBackedLikFitMixin`

   ### 2. `InferenceAbstractKKGEE`

   This class is also a good MixIn candidate.

   It packages together:

   - KK design eligibility checks
   - GEE model generation
   - extraction of coefficient / variance structures
   - GEE-specific helper methods and fallback logic

   Conceptually this is closer to a `KKGEEEngineMixin` than to a true family in the semantic inheritance tree.

   ### 3. `InferenceAbstractKKGLMM`

   Same reasoning as GEE:

   - KK design restrictions
   - GLMM backend setup
   - family-specific abstract hooks such as `glmm_response_type()` and `glmm_family()`
   - shared extraction / fit logic

   This looks more like a reusable engine bundle than a semantic superclass.

   ### 4. `InferenceKKPassThrough`

   This is a borderline MixIn candidate.

   It currently bundles:

   - KK design validation
   - a “pass through to downstream model/test implementation” style
   - one-likelihood / IVWC helper infrastructure

   The reason it is only borderline is that the package currently treats `InferenceKKPassThrough` as a meaningful category in several places. If it became a MixIn, we would need an explicit replacement for that capability detection.

### 5. `InferenceKKPassThroughCompound`

This is another good behavioral candidate:

- compound KK passthrough logic
- shared helpers for robust regression, quantile, OLS-one-likelihood, Bai-adjusted variants

It is a feature bundle more than a semantic base.

### 6. `InferenceAbstractQuantileRandCI`

This is a good narrow MixIn candidate for:

- quantile-specific randomization-CI behavior
- the Zhang-combined-CI machinery constraints

This is very behavior-specific and not an especially important semantic node.

## Classes that are only weak MixIn candidates

These could be partially decomposed, but I would not convert the entire class into a MixIn immediately.

### `InferenceParamBootstrap`

Right now it is thin enough to be either:

- a true inheritance node, or
- a future capability marker / MixIn

But because it is intended to mean something package-wide, keeping it as a real superclass is cleaner unless the whole inference hierarchy is refactored around explicit capability roles.

### `InferenceMLEorKMSummaryTable`

This is more of a specialized concrete-ish infrastructure class than a broad reusable bundle. It is not the first place I would introduce MixIns.

## What would likely become MixIns in a refactor

If we moved to MixIn-style composition, the most likely conversions are:

- `InferenceAsympLikStdModCache`
- `InferenceAbstractKKGEE`
- `InferenceAbstractKKGLMM`
- `InferenceKKPassThrough`
- `InferenceKKPassThroughCompound`
- `InferenceAbstractQuantileRandCI`

Potentially also sub-bundles extracted from those classes, which would be even cleaner:

- `KKDesignEligibilityMixin`
- `StdModelCacheMixin`
- `LikelihoodTestSpecDispatchMixin`
- `DesignBackedBootstrapWorkerMixin`
- `GEEEngineMixin`
- `GLMMEngineMixin`
- `QuantileRandCIMixin`

That smaller-grained approach is probably better than converting the current large abstract classes wholesale into MixIns.

## Would MixIns allow the current `is(...)` checks?

Short answer: **no, not automatically**.

If a MixIn is implemented by copying methods into a class definition, then:

- `inherits(self, "SomeMixin")` will be `FALSE`
- `is(self, "SomeMixin")` will not become `TRUE`
- `generator$get_inherit()` will not show the MixIn
- `InferenceSuite`-style inheritance walking will not see the MixIn

That is because the MixIn is not part of the real R6 superclass chain.

## What would still work

These would still work:

- `class(self)[1]` for the concrete class name
- direct method calls added by the MixIn
- explicit marker fields or capability methods we define ourselves

## What would have to change

If we adopt MixIns, the package should stop using superclass membership as a proxy for capabilities.

Instead of:

- `inherits(self, "InferenceAsympLikStdModCache")`
- `inherits(self, "InferenceAbstractKKGEE")`
- `inherits(self, "InferenceKKPassThrough")`

we would want explicit capability checks such as:

- `private$has_role("std_model_cache")`
- `private$has_role("gee_engine")`
- `private$has_role("kk_passthrough")`
- `self$supports_likelihood_tests()`
- `self$supports_quantile_rand_ci()`

or a static class-level marker like:

- `private$capabilities = c("std_model_cache", "design_backed_bootstrap", "gee_engine")`

## Repo-specific impact on current checks

### `InferenceSuite`

Current logic in [EDI/R/inference_suite.R](../EDI/R/inference_suite.R) walks `get_inherit()`.

That means:

- MixIns would not be discovered as ancestors
- “abstract base class” exclusion by class name would no longer be enough

This file would need a redesign if we want suite discovery to understand MixIn-composed infrastructure classes.

### `InferenceRandCI`

Current logic in [EDI/R/inference_all_abstract_rand_ci.R](../EDI/R/inference_all_abstract_rand_ci.R) uses a long `inherits(...)` list to determine whether to apply GLM-like transforms in randomization CI calculations.

This is exactly the kind of logic that MixIns make worse unless replaced by explicit capability predicates.

### Other abstract-family branch checks

The quantile branch checks in [EDI/R/inference_all_abstract_rand.R](../EDI/R/inference_all_abstract_rand.R) would also need to become explicit capability checks rather than superclass tests.

## Recommended design if we pursue MixIns

I would not do a full “replace abstract inheritance with MixIns” rewrite.

I would do this instead:

1. Keep the semantic superclass spine intact:
   - `Inference`
   - `InferenceRand`
   - `InferenceRandCI`
   - `InferenceNonParamBootstrap`
   - `InferenceAsymp`
   - `InferenceAsympLik`
   - `InferenceExact`

2. Introduce small MixIn helpers only for reusable implementation bundles:
   - standardized model-cache behavior
   - design-backed bootstrap-worker behavior
   - GEE engine behavior
   - GLMM engine behavior
   - KK design eligibility checks
   - quantile randomization-CI behavior

3. Add an explicit capability system:
   - `has_role()`
   - `supports_*()` methods
   - or static role vectors on each class

4. Gradually replace `inherits(self, "...")` checks that are really capability checks.

## Bottom line

MixIns would help most with code reuse for classes like:

- `InferenceAsympLikStdModCache`
- `InferenceAbstractKKGEE`
- `InferenceAbstractKKGLMM`
- `InferenceKKPassThrough`
- `InferenceKKPassThroughCompound`
- `InferenceAbstractQuantileRandCI`

But MixIns would **not** preserve the current `is(...)`, `inherits(...)`, or `get_inherit()` behavior by themselves.

So the package can benefit from MixIns only if we also stop treating inheritance as the main capability-discovery mechanism.

Without that additional redesign, strict inheritance remains the safer architecture for this repo.
