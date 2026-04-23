# Information Preference Defaults

The public default is always `information_preference = "auto"`.

`"auto"` resolves through the internal `get_default_information_source()` hook:

1. `"fisher"` when the inference path exposes a genuine analytic Fisher information matrix.
2. `"observed"` when the path exposes only observed information.
3. `"legacy"` only as a fallback for older information-backed paths that still expose only a generic `information` matrix.

The most recent source actually used is available through `get_information_source_used()`.

## Current `auto -> fisher` paths

These paths currently advertise `supports_fisher_information() = TRUE`:

- `InferenceIncidLogRegr`
- `InferenceCountPoisson`
- `InferenceIncidKKClogitOneLik`
- `InferenceCountPoissonUnivKKCPoissonCombinedLikelihood`
- `InferenceCountPoissonMultiKKCPoissonCombinedLikelihood`

## Current `auto -> observed` paths

These families currently support information-backed inference but do not expose analytic Fisher information:

- KK GLMM logistic/proportion families via `InferenceAbstractKKGLMM`
- KK clogit-plus-GLMM one-likelihood families via `InferenceAbstractKKClogitPlusGLMM`
- zero-augmented / hurdle count families via `InferenceCountZeroAugmentedPoisson`
- KK hurdle-Poisson one-likelihood families via `InferenceAbstractKKHurdlePoissonOneLik`
- KK Weibull frailty one-likelihood families via `InferenceAbstractKKWeibullFrailtyOneLik`

## Notes

- `wald`, `score`, and future information-based tests should all consult the same information source selection.
- Wald-only classes that do not expose switchable information-backed inference still support only `"auto"`.
- When a class exposes both observed and Fisher information, `"auto"` should stay stable unless the class-specific default hook changes.
