# InferenceAsympLik Design Note

## Proposal

Add an intermediate `InferenceAsympLik` class between `InferenceAsymp` and the likelihood-backed families:

- `InferenceMLEorKMforGLMs`
- `InferenceKKPassThrough`
- `InferenceAbstractKKGLMM`

but not the other asymptotic families such as bootstrap, exact, randomization-only, GEE, or g-computation classes.

## Short Answer

This is a reasonable refactor, but I would treat it as a cleanup pass, not as a correctness fix.

The current code already has a working likelihood-backed abstraction through `InferenceAsymp`, and the three candidate branches already use the same asymptotic test machinery:

- Wald p-values and intervals
- score p-values and intervals
- gradient p-values and intervals
- likelihood-ratio p-values and intervals
- likelihood-test memoization and warm-start plumbing

So the main value of `InferenceAsympLik` would be organizational:

- move likelihood-specific helper methods out of `InferenceAsymp`
- keep `InferenceAsymp` thinner
- make the inheritance tree more explicit

That is useful, but only if the move is done carefully.

## Why It Makes Sense

There is a real conceptual split in the current hierarchy:

- `InferenceAsymp` is used by many families that only need Wald-style asymptotics.
- a smaller subset also exposes a likelihood specification and can support score / gradient / likelihood-ratio paths.

The likelihood-backed branches all need the same shared scaffolding:

- `get_likelihood_test_spec()`
- `get_standard_error()`
- `get_degrees_of_freedom()`
- `compute_score_two_sided_pval_impl()`
- `compute_gradient_two_sided_pval_impl()`
- `compute_lik_ratio_two_sided_pval_impl()`
- likelihood-test caching and warm starts
- score/information extraction
- CI inversion for score, gradient, and likelihood ratio

That is enough commonality to justify an intermediate base class.

## Why I Would Be Cautious

`InferenceAsymp` is already doing more than “plain asymptotic” work. It currently mixes:

- Wald dispatch
- likelihood-test dispatch
- memoization
- CI inversion
- information-matrix selection
- warm-start state

So introducing `InferenceAsympLik` only helps if you are willing to move the likelihood-specific block out of `InferenceAsymp` cleanly. If you do not, you just add another layer without reducing complexity.

There is also a maintenance risk:

- method lookup in R6 is sensitive to where the concrete override lives
- some subclasses already override `compute_asymp_confidence_interval()` or `shared()`
- the KK GLMM classes and the MLE/KM classes use different fitter conventions

A shallow split can become a source of confusion if the base classes overlap too much.

## My Recommendation

I would do it only if the goal is broader architectural cleanup, not just to introduce a new name.

Best version:

1. Keep `InferenceAsymp` as the generic Wald/asymptotic shell.
2. Move the likelihood-specific machinery into `InferenceAsympLik`.
3. Make `InferenceMLEorKMforGLMs`, `InferenceKKPassThrough`, and `InferenceAbstractKKGLMM` inherit from `InferenceAsympLik`.
4. Leave non-likelihood families on `InferenceAsymp`.

That would make the class tree clearer and reduce the impression that every asymptotic model must know about likelihood inversion internals.

## Cost / Risk Estimate

This is a medium refactor, not a trivial rename.

Main risks:

- regression in score / gradient / likelihood-ratio CI inversion
- accidental breakage of classes that inherit `InferenceAsymp` but do not currently advertise likelihood tests
- duplicated override logic in KK special cases
- test churn across many families

Expected work:

- low if this is just a structural extraction with no semantic change
- moderate if you want the class split to be clean and complete
- high if you also want to simplify the subclass overrides at the same time

## Bottom Line

I think `InferenceAsympLik` is a good idea if you want the inheritance tree to reflect the distinction between:

- “generic asymptotic inference”
- “asymptotic inference with likelihood-based testing”

But I would not do it just to shorten `InferenceAsymp`. The refactor only pays off if you also use it to isolate the likelihood-testing machinery that is now spread across the current base and the three main descendant branches.
