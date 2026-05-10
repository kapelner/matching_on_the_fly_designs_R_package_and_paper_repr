# Count Likelihood Inheritance Proposal

## Goal

Split the count-model likelihood code away from `InferenceMLEorKMforGLMs` without losing the shared likelihood-test and fit-reuse machinery.

The current situation mixes two different concepts:

- standard GLM / KM likelihood-backed inference
- count-specific likelihood models with custom parameterizations

Those are similar enough to share the asymptotic test machinery, but they are not the same abstraction.

## Proposed Hierarchy

### 1. Keep the existing GLM/KM branch

No structural change:

```text
InferenceAsympLik
└─ InferenceMLEorKMforGLMs
   ├─ InferenceCountPoisson
   ├─ InferenceCountNegBin
   ├─ InferenceIncidLogRegr
   ├─ InferenceIncidLogBinomial
   ├─ InferenceIncidModifiedPoisson
   ├─ InferencePropBetaRegr
   ├─ InferenceSurvivalCoxPHRegr
   ├─ InferenceSurvivalWeibullRegr
   └─ many other standard GLM / KM classes
```

This class should remain the shared base for the ordinary regression family that already fits the `generate_mod() -> beta_hat_T / s_beta_hat_T / df` pattern.

### 2. Add a count-likelihood base

New intermediate class:

```text
InferenceAsympLik
└─ InferenceCountLikelihood
   ├─ InferenceCountZeroAugmentedPoissonAbstract
   │  ├─ InferenceCountZeroInflatedPoisson
   │  ├─ InferenceCountZeroInflatedNegBin
   │  └─ InferenceCountHurdlePoisson
   ├─ InferenceCountHurdleNegBin
   └─ InferenceCountPoisson / InferenceCountNegBin
```

This new base should own the count-specific shared logic that is currently duplicated or awkwardly routed through the GLM/KM base:

- count-specific parameter packing and unpacking
- count-specific warm starts
- count-specific null-fit construction
- score / information / loglik dispatch for count models
- estimate, SE, and df caching for count likelihood families

### 3. Add a companion-model branch for robust / quasi count models

Separate new base:

```text
InferenceAsympLik
└─ InferenceCountCompositeLikelihood
   ├─ InferenceCountRobustPoisson
   └─ InferenceCountQuasiPoisson
```

This should not pretend to be a full likelihood base for the reported estimator. It should expose a companion Poisson likelihood solely for test inversion and p-value generation, while keeping the reported inference robust or quasi.

## Exact Parent / Child Moves

### Move to `InferenceCountLikelihood`

| Current class | Proposed parent |
|---|---|
| `InferenceCountPoisson` | `InferenceCountLikelihood` |
| `InferenceCountNegBin` | `InferenceCountLikelihood` |
| `InferenceCountZeroAugmentedPoissonAbstract` | `InferenceCountLikelihood` |
| `InferenceCountHurdleNegBin` | `InferenceCountLikelihood` |

### Keep under `InferenceCountZeroAugmentedPoissonAbstract`

| Current class | Proposed parent |
|---|---|
| `InferenceCountZeroInflatedPoisson` | `InferenceCountZeroAugmentedPoissonAbstract` |
| `InferenceCountZeroInflatedNegBin` | `InferenceCountZeroAugmentedPoissonAbstract` |
| `InferenceCountHurdlePoisson` | `InferenceCountZeroAugmentedPoissonAbstract` |

### Move to `InferenceCountCompositeLikelihood`

| Current class | Proposed parent |
|---|---|
| `InferenceCountRobustPoisson` | `InferenceCountCompositeLikelihood` |
| `InferenceCountQuasiPoisson` | `InferenceCountCompositeLikelihood` |

## Why This Split Is Better

`InferenceMLEorKMforGLMs` already does useful work and should stay in place. The count classes are different enough that forcing them through that class makes the abstraction leaky:

- zero-augmented models need two linear predictors and special null-fit handling
- negative binomial models carry an extra dispersion parameter
- robust/quasi count classes are not full likelihood estimators at all, but still want a likelihood-analogue for tests

That means a count-specific layer is the right place to centralize the count mechanics without contaminating the GLM/KM layer.

## Migration Order

### Phase 1: Introduce `InferenceCountLikelihood`

Start with the models that are already closest to full-likelihood count regression:

1. `InferenceCountPoisson`
2. `InferenceCountNegBin`

These two are the cleanest proof that the new base is viable.

### Phase 2: Move zero-augmented count models

3. `InferenceCountZeroAugmentedPoissonAbstract`
4. `InferenceCountZeroInflatedPoisson`
5. `InferenceCountZeroInflatedNegBin`
6. `InferenceCountHurdlePoisson`

This is the point where the separate count-specific null-fit and warm-start logic becomes obviously worth the new layer.

### Phase 3: Move the remaining true count-likelihood model

7. `InferenceCountHurdleNegBin`

This class has bespoke likelihood machinery, but it still belongs in the count-likelihood family rather than the GLM/KM family.

### Phase 4: Add the companion-model branch

8. `InferenceCountRobustPoisson`
9. `InferenceCountQuasiPoisson`

Do these last. They are structurally different because the likelihood path is a companion model, not the reported estimator.

## Recommendation

I would not try to flatten `InferenceMLEorKMforGLMs` into the count branch.

Instead:

- keep `InferenceMLEorKMforGLMs` for ordinary GLM/KM models
- add `InferenceCountLikelihood` for true count-likelihood models
- add `InferenceCountCompositeLikelihood` for robust/quasi companion-model testing

That gives you a cleaner hierarchy without forcing incompatible model families through the same abstraction.
