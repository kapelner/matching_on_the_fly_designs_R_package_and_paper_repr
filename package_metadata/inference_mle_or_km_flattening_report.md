# InferenceMLEorKMforGLMs Flattening Report

## Conclusion

`InferenceMLEorKMforGLMs` is not duplicative in the sense of a removable wrapper. It is a real shared implementation base for the package’s likelihood-backed GLM and KM-style inference classes. I would **not** flatten it.

The class owns shared fit caching, estimate extraction, bootstrap worker reuse, randomization-inference reuse, standard-error extraction, degrees-of-freedom handling, and likelihood-test dispatch. If it were flattened into each concrete subclass, that logic would have to be copied into many children or moved into yet another replacement base class with effectively the same responsibility.

## What It Provides

| Responsibility | Present in `InferenceMLEorKMforGLMs` | Why it matters |
|---|---:|---|
| Shared treatment-estimate extraction | Yes | Concrete subclasses rely on `shared()` to populate `beta_hat_T`. |
| Shared bootstrap worker path | Yes | Many subclasses reuse the design-backed bootstrap state. |
| Shared randomization-inference path | Yes | `compute_treatment_estimate_during_randomization_inference()` is centralized here. |
| Standard-error extraction | Yes | The class chooses between information-based and model-based SE paths. |
| Degrees of freedom caching | Yes | Several concrete classes use t-calibrated inference through this layer. |
| Likelihood-test plumbing | Yes | Score, gradient, and likelihood-ratio paths are dispatched here. |
| Warm-start and null-fit memoization | Yes | The likelihood inversion logic depends on this shared machinery. |

## Evidence It Is Not Just a Wrapper

The source file [EDI/R/inference_all_abstract_mle_or_KM_for_GLMs.R](../EDI/R/inference_all_abstract_mle_or_KM_for_GLMs.R) contains substantial private implementation:

- `shared()`
- `generate_mod()` as an abstract hook used by subclasses
- `create_bootstrap_worker_state()`
- `load_bootstrap_sample_into_worker()`
- `compute_bootstrap_worker_estimate()`
- `get_standard_error()`
- `get_degrees_of_freedom()`
- `compute_score_two_sided_pval_impl()`
- `compute_gradient_two_sided_pval_impl()`
- `compute_lik_ratio_two_sided_pval_impl()`
- `make_warm_fit_null_wrapper()`
- `compute_likelihood_test_two_sided_pval()`

That is not inheritance noise. It is the core of the model-fitting and asymptotic test infrastructure for the GLM/KM branch.

## Current Children

Representative direct children include:

- `InferenceCountPoisson`
- `InferenceCountNegBin`
- `InferenceCountRobustPoisson`
- `InferenceCountQuasiPoisson`
- `InferenceCountZeroInflatedPoisson`
- `InferenceCountZeroInflatedNegBin`
- `InferenceCountHurdlePoisson`
- `InferenceCountHurdleNegBin`
- `InferenceIncidLogBinomial`
- `InferenceIncidModifiedPoisson`
- `InferencePropBetaRegr`
- `InferencePropFractionalLogit`
- `InferencePropZeroOneInflatedBetaRegr`
- `InferenceSurvivalCoxPHRegr`
- `InferenceSurvivalWeibullRegr`
- `InferenceSurvivalDepCensTransformRegr`
- `InferenceOrdinalAdjCatLogitRegr`
- `InferenceOrdinalCauchitRegr`
- `InferenceOrdinalCloglogRegr`
- `InferenceOrdinalOrderedProbitRegr`
- `InferenceOrdinalPropOddsRegr`

That breadth is the point of the class: it prevents the same fit/SE/test code from being duplicated across many model families.

## What Could Be Flattened Instead

There is one nearby abstraction that still looks flattenable:

| Parent | Child | Status |
|---|---|---|
| `InferenceSurvivalStratCoxPHAbstract` | `InferenceSurvivalStratCoxPHRegr` | Candidate for flattening |

That is a separate cleanup from `InferenceMLEorKMforGLMs`. The survival stratified Cox abstract is a single-child wrapper on top of this base class.

## Recommendation

Keep `InferenceMLEorKMforGLMs` as a shared implementation base.

If the goal is hierarchy cleanup, focus on the single-child abstract wrappers above it, not this class. Flattening `InferenceMLEorKMforGLMs` would remove a useful seam and replace it with duplicated code or a new equivalent base under another name.
