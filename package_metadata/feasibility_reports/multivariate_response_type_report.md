# Multivariate / Vector-Valued Response Type Report

## Scope

This report evaluates how difficult it would be to add first-class support to
the package for **multivariate / vector-valued outcomes**: experiments where
each subject contributes more than one outcome measurement that should be
analyzed jointly rather than one-at-a-time, for example:

- co-primary clinical endpoints
- a biomarker panel measured once per subject
- a bundle of educational or social outcomes
- a joint physics or chemistry response surface

This report was commissioned by
[response_types_landscape_report.md](references/response_types_landscape_report.md)'s
TODO-4, which flagged this family as `very hard` based on the "Clinical
trials" field section's co-primary-endpoint and biomarker-panel examples. That
document only motivates the gap qualitatively; this report does the concrete
architectural assessment.

This is a fundamentally different kind of extension than adding one more
`response_type` value (e.g. the `nominal` case covered in
[nominal_response_type_report.md](nominal_response_type_report.md)). Adding
`nominal` still fits one scalar response per subject into the existing
per-subject storage; multivariate does not.

## Short Answer

The verdict depends entirely on which of two very different features is
meant by "multivariate support," and the landscape survey's single `very
hard` label does not distinguish them:

1. **Composite/multiplicity-adjusted analysis of several existing univariate
   outcomes fit independently, then combined** — **moderate**. This requires
   no change to the one-scalar-response-per-subject core; it is a new
   orchestration layer on top of what already exists.
2. **True joint/multivariate modeling** (a single model that estimates a
   vector-valued treatment effect with a joint covariance, e.g. seemingly
   unrelated regression, multivariate GLM, or a global rank-sum-type test on
   a response vector) — **very hard**, and correctly so. This requires new
   per-subject response storage, new C++ numerical cores, and a redesign of
   every layer that currently assumes one scalar estimand.

So `very hard` is the right verdict for genuine joint modeling, but it
overstates the cost of the composite/multiplicity path, which is the
practical way most applied co-primary-endpoint analyses are actually done
(run each endpoint's own model, then control the family-wise error rate
across the K tests).

## How Common Are Multivariate/Co-Primary Outcomes In Experimental Literatures?

- **Clinical trials**: co-primary endpoints and biomarker panels are common
  enough that ICH E9 and FDA multiplicity guidance both address them
  explicitly; regulatory trials frequently pre-specify 2-4 co-primary or key
  secondary endpoints with a formal multiplicity-control plan (e.g.
  hierarchical testing, Bonferroni/Holm, or a global test) rather than a true
  joint parametric model.
- **Education and social science**: outcome "bundles" (e.g. multiple test
  subscales, or several behavioral measures) are common, and applied
  practice there also leans heavily on per-outcome testing with a
  multiple-comparisons correction, or an index/composite score computed
  before the analysis (which then re-enters the package as an ordinary
  scalar `continuous` outcome — no new machinery needed).
- **True joint parametric multivariate modeling** (vector-valued treatment
  effects with an explicit cross-endpoint covariance, e.g. SUR-style
  regression) is comparatively rare in the applied RCT/field-experiment
  literature relative to the composite/multiplicity-adjusted approach, largely
  *because* it is harder to interpret, harder to communicate to
  non-statistician stakeholders, and requires larger samples to estimate the
  cross-endpoint covariance well.

Practical conclusion: the dominant applied need is composite/multiplicity
analysis of several already-supported scalar outcomes, not full joint
modeling. This should reorder the package's priorities relative to a naive
reading of the landscape survey's "very hard" tag.

## What Already Exists

Confirmed via grep of `EDI/R/*.R`:

- There is **no multiplicity-correction infrastructure anywhere in the
  package** (`grep -rn "multiplicity\|multiple_testing\|holm\|bonferroni\|p\.adjust" EDI/R/*.R`
  returns nothing except one unrelated docstring mention of
  "multiplicity-corrected" in
  [simulation_framework_report.R:83](../EDI/R/simulation_framework_report.R:83)).
  So even the "easier" composite path (#1 above) requires new code, not
  reuse of an existing correction layer — it is just much less new code than
  #2.
- There is **no multi-column response storage** anywhere in `Design`. The
  response is stored as a single numeric vector `private$y`, written by
  `add_one_subject_response(t, y, dead)` (`y` is asserted `len = 1`
  numeric — [design_abstract.R:145-148](../EDI/R/design_abstract.R:145)) and
  `add_all_subject_responses(ys, deads)` (`private$y = as.numeric(ys)` —
  [design_abstract.R:200-229](../EDI/R/design_abstract.R:200)). This is the
  central architectural fact this report is built around: every experiment
  currently has exactly one `Design` object holding exactly one outcome
  vector.
- `response_type` validation is a fixed `assertChoice()` over six values —
  `continuous`, `incidence`, `proportion`, `count`, `survival`, `ordinal`
  ([design_abstract.R:94](../EDI/R/design_abstract.R:94)) — with no vector or
  matrix-valued option, and no `nominal` either (confirming that report's own
  premise that nominal hasn't landed yet).
- `SimulationFramework`'s `betaT` parameter accepts "a numeric scalar or
  vector" ([simulations_framework.R:329](../EDI/R/simulations_framework.R:329)),
  but this is a **sweep of scalar effect-size values across separate DGP
  cells** (`betaT_values = unique(as.numeric(betaT))` at
  [simulations_framework.R:625](../EDI/R/simulations_framework.R:625)), not a
  single vector-valued treatment effect applied jointly to several outcomes
  in one replication. Each cell still simulates and fits one scalar outcome.
  This is worth flagging explicitly because it is easy to misread as
  existing multivariate support — it is not.

## Difficulty By Layer

### 1. Design base class: hard for joint modeling, unnecessary for composite path

For the composite path (#1), no `Design`-layer change is needed at all: each
of the K endpoints is simply its own `Design` object (or the same covariates
threaded through K parallel `Design` objects, one per outcome), fit with the
existing single-response machinery K times.

For true joint modeling (#2), `Design` would need:
- storage for a response **matrix** (`n × K`) instead of a response vector
- a `response_type` per column (endpoints can have different families —
  e.g. one continuous biomarker and one binary responder-status endpoint)
- assertion/validation logic for partial missingness per endpoint (a subject
  present in the trial but missing one of several co-primary measurements)

This is a genuine breaking change to the `Design` contract, not an additive
field, since `add_one_subject_response`/`add_all_subject_responses` are typed
around a single scalar `y` throughout.

### 2. Most design classes: unaffected for composite path, hard for joint path

Design classes that are response-type agnostic (assignments, covariates,
timing) are untouched by the composite path — they are simply instantiated
K times, once per endpoint. For the joint path, any design class with
response-adaptive logic (see next point) would need matrix-aware treatment.

### 3. KK21 / response-adaptive weighting designs: moderate for composite (reuse per-endpoint), hard for joint

`DesignSeqOneByOneKK21` computes covariate weights from the observed response
family with explicit per-`response_type` branches
([design_seq_one_by_one_KK21.R:253-281](../EDI/R/design_seq_one_by_one_KK21.R:253)).
For composite analysis this is unaffected (each endpoint's design/weighting
is independent). For joint modeling, response-adaptive weighting would need
to either pick one "primary" endpoint to drive allocation (a reasonable and
common simplification) or genuinely fuse multiple endpoints into one weight
— the latter is a substantial, unresolved design question, not just an
implementation task.

### 4. `SimulationFramework`: moderate for composite, hard for joint

`SimulationFramework` currently transforms one latent continuous signal
`y_cont` to one response-family scale via a single
`switch(response_type, ...)` dispatcher
([simulations_framework.R:126](../EDI/R/simulations_framework.R:126)) and
defines response-type-specific treatment-effect semantics and curated
default inference classes per response type
([simulations_framework.R:3106](../EDI/R/simulations_framework.R:3106)).

For the composite path, `SimulationFramework` could run K independent
single-outcome simulations (potentially with correlated latent draws across
endpoints, which is a moderate addition: draw a correlated multivariate
normal latent vector, then apply the existing per-endpoint
`transform_cont_y_based_on_response_type()` column-wise — this function
already exists standalone at
[simulations_framework.R:105-124](../EDI/R/simulations_framework.R:105) and
is reusable as-is per endpoint).

For joint modeling, `SimulationFramework` would need a vector-valued `betaT`
applied *within* one replication (not swept across cells), a
vector/matrix-valued truth object, and vector-valued MSE/coverage/power
summaries — none of which the current scalar-truth summary machinery
supports.

## The Main Conceptual Problem: What Is The Multivariate Estimand?

Exactly as with nominal outcomes, this is the crux. Unlike nominal (where the
ambiguity is about *which contrast*), here the ambiguity is about *whether
there should be one number at all*:

- **One global test statistic** (e.g. Hotelling-type or O'Brien-type combined
  test): answers "is there evidence of a treatment effect somewhere across
  the endpoint panel," but is not, by itself, an effect estimate with a CI
  for any specific quantity.
- **K per-endpoint scalar estimates with a multiplicity-adjusted decision
  rule**: gives K familiar scalar estimates (reusing 100% of the existing
  per-endpoint machinery), with a shared correction layer (Holm/Bonferroni,
  or a hierarchical/gatekeeping procedure) deciding which are declared
  significant. This is by far the best fit for the package's current
  one-scalar-estimand architecture.
- **A single joint vector-valued treatment effect with a joint covariance
  matrix**: statistically the "purest" version, but breaks essentially every
  scalar-estimand assumption described below.

## How The Existing `inference_all_*` Paths Would Respond

### The scalar-treatment-effect assumption is pervasive

`InferenceMLEorKMSummaryTable` and `InferenceAsympLikStdModCache`, the shared
abstract layers underneath most likelihood-based inference paths, both
assume a single scalar treatment coefficient (`generate_mod()` returning an
object with `b[2]` as the treatment effect —
[inference_all_abstract_asymp_lik_std_mod_cache.R:6](../EDI/R/inference_all_abstract_asymp_lik_std_mod_cache.R:6)),
and the base caching layer in `inference_all_abstract.R` stores exactly one
`beta_hat_T` / `s_beta_hat_T` pair per fitted object
([inference_all_abstract.R:592-602](../EDI/R/inference_all_abstract.R:592)).
None of this can hold a vector-valued treatment effect without a structural
change to the cache contract that essentially every concrete `Inference*`
class relies on.

### The composite path needs none of that changed

Because the composite path fits K independent scalar models (one call into
the existing `Inference*` machinery per endpoint), it needs zero changes to
`inference_all_abstract.R`, `InferenceAsympLikStdModCache`, or any concrete
inference class. The only new code is a thin orchestration wrapper that:
1. runs the same (or different) `Inference*` class once per endpoint column,
2. collects the K `(beta_hat_T, s_beta_hat_T, pval)` triples,
3. applies a multiplicity correction across the K p-values,
4. optionally computes a combined global-test statistic from the K results.

### `InferenceSuite`: works for the composite path, not applicable to joint modeling

`InferenceSuite` discovers applicable classes by trying to instantiate each
one against a `Design` and reading off compatibility
([inference_suite.R:110](../EDI/R/inference_suite.R:110)). For the composite
path, `InferenceSuite` is unaffected: it would simply be invoked K times,
once per per-endpoint `Design`. For joint modeling, `InferenceSuite`'s
per-class discovery model would need to understand multi-endpoint `Design`
objects, which it currently has no concept of.

## Package Areas That Need Explicit Guards (Joint-Modeling Path Only)

If a matrix-valued `Design` were ever introduced for true joint modeling, the
following would need explicit new-response-shape assertions, mirroring the
nominal report's "fence off generic paths" recommendation:
- `InferenceAllSimpleMeanDiff` / `InferenceAllSimpleMeanDiffPooledVar` /
  `InferenceAllSimpleWilcox` — all currently assume `private$y` is a plain
  numeric vector; a matrix response would silently break their internal
  arithmetic (mean of a matrix column vs. mean of the whole matrix) rather
  than erroring cleanly, which is worse than nominal's failure mode (nominal
  at least produces a runnable-but-meaningless numeric answer; a matrix `y`
  passed into code written for a vector could error unpredictably deep
  inside `mean()`/`var()` calls, or silently vectorize incorrectly).
- `Design$assert_y()` ([design_abstract.R:495](../EDI/R/design_abstract.R:495))
  would need a whole new branch for matrix-shaped `y`.

The composite path needs no guards anywhere, since it never constructs a
matrix-valued `Design`.

## Recommended First-Wave Multivariate Support

Consistent with the "Short Answer" above, the pragmatic first wave is
entirely the composite path:

1. `InferenceMultiEndpointComposite` — a thin R6 orchestration class that
   takes K `(Design, Inference class, ctor args)` tuples sharing the same
   assignment vector `w`, runs each independently, and returns a summary
   table of K estimates plus multiplicity-adjusted p-values (Holm as the
   default, since it strictly dominates Bonferroni and needs no
   distributional assumption beyond what each endpoint's own p-value already
   has).
2. A global combined test (e.g. Fisher's combination of the K per-endpoint
   p-values, or an O'Brien-type combined-rank test if the endpoints share a
   common direction-of-benefit convention) as a second, optional output of
   the same class — this gives the "is there a signal anywhere" answer
   without requiring a joint covariance estimate.
3. Defer true joint/vector-valued modeling (SUR-style regression, a
   Hotelling-type joint Wald test, or joint bootstrap over a response
   matrix) to a second-generation project, gated on whether a concrete user
   need for a *joint* covariance estimate (as opposed to per-endpoint
   estimates plus a combined test) actually arises. Most applied
   co-primary-endpoint practice does not require it.

## `SimulationFramework` Implications

### Composite path

Feasible without breaking any existing scalar-truth contract:
- generate K correlated latent signals (a multivariate normal draw with a
  user-specified cross-endpoint correlation matrix, reusing
  `transform_cont_y_based_on_response_type()` column-wise, per endpoint),
- run K independent fits per replication via the composite wrapper above,
- report K sets of the existing scalar MSE/coverage/power summaries, plus
  one new family-wise power/type-I-error summary across the K decisions
  (whether the multiplicity-adjusted procedure correctly detects/rejects).

### Joint-modeling path

Would require a genuinely new truth object (`betaT` as a vector, applied
within one replication), a new vector/matrix-valued MSE-analogue (e.g.
Mahalanobis-distance-based), and new curated default inference classes keyed
on "multivariate" the same way
[simulations_framework.R:3106](../EDI/R/simulations_framework.R:3106) keys
scalar defaults on `response_type` today. This is the single most expensive
integration point for the joint path, exactly as `SimulationFramework` was
the hardest integration point for nominal.

## Recommended Implementation Plan

### Stage 1: Composite orchestration layer

Add `InferenceMultiEndpointComposite` as a wrapper over existing single-
endpoint `Inference*` classes, with a Holm-adjusted p-value table. Requires
no `Design` changes.

### Stage 2: Combined global test

Add one combined-test option (Fisher combination or O'Brien-type rank
combination) to the same wrapper, giving a scalar "any endpoint affected"
answer without a joint covariance estimate.

### Stage 3: `SimulationFramework` composite support

Extend `SimulationFramework` to generate K correlated latent outcomes and
drive the Stage 1/2 wrapper, with family-wise power/type-I-error summaries.

### Stage 4 (only if a concrete need emerges): true joint modeling

Introduce a matrix-valued `Design` variant, new C++ numerical cores for
joint/SUR-style estimation, and the accompanying `SimulationFramework`
vector-truth machinery. This stage should not be started speculatively —
gate it on an actual user request for a joint covariance estimate, since the
composite path already serves the dominant applied use case at a fraction of
the cost.

## Bottom Line

The landscape survey's single `very hard` verdict is correct only for **true
joint/vector-valued modeling**, which does require breaking the package's
one-scalar-response-per-subject `Design` contract, new C++ numerical cores,
and a redesign of `SimulationFramework`'s truth/summary machinery — a
second-generation architectural project on par with (or larger than) the
`SimulationFramework` cost identified in the nominal-response-type report.

But the dominant applied version of "multivariate support" — composite
analysis of several co-primary endpoints with a multiplicity-adjusted
decision rule, which is how most real clinical and field-experiment
co-primary-endpoint analyses are actually run — is **moderate**, requires
zero changes to `Design`, `inference_all_abstract.R`, or any concrete
`Inference*` class, and can be built entirely as a new orchestration layer
plus one new multiplicity-correction utility (which does not exist anywhere
in the package today).

Recommendation: implement the composite path first (Stages 1-3). Do not
commission a `Design`-matrix / joint-modeling project (Stage 4) until a
concrete use case demands a genuine joint covariance estimate rather than
per-endpoint estimates plus a combined decision rule.
