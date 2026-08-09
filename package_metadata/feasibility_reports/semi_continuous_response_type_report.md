# Semi-Continuous / Zero-Inflated Continuous Response Type Report

## Scope

This report evaluates how difficult it would be to add first-class support for
**semi-continuous outcomes** to the package: a continuous response with a
point mass at a boundary (almost always zero) plus a continuous distribution
elsewhere. Typical examples:

- medical costs / health-care utilization expenditure (many subjects spend $0)
- environmental or laboratory measurements with a detection floor
- treatment intensity or dosage with many true zeros (never treated)

This fulfills TODO-7 in
[references/response_types_landscape_report.md](references/response_types_landscape_report.md),
which flagged this family as `moderate to hard` and specifically noted that
EDI already has zero-inflated/two-part paths for two *other* response
families — `proportion` (`InferencePropZeroOneInflatedBetaRegr`) and `count`
(the ZIP/ZINB/Hurdle-Poisson/Hurdle-NegBin family) — but none for
`continuous`.

The report covers two distinct sub-cases that the landscape survey collapsed
into one difficulty label, because they turn out to warrant different
verdicts:

1. **Structural zero-inflation** (Tobit-like): the response is exactly zero
   for a subset of subjects for a real, modelable reason (no cost incurred,
   no dose given), and continuously distributed and positive otherwise.
2. **Detection-limit censoring**: the response is *unobserved* below some
   threshold, not truly zero — a censoring problem, not a mixture problem.

## Short Answer

Adding a **structural semi-continuous / hurdle-continuous inference path**
that reuses the existing count zero-augmented architecture is **easy to
moderate** — the single easiest of the six architectural-extension families
evaluated in this batch of reports (alongside
[longitudinal_repeated_measures_response_type_report.md](longitudinal_repeated_measures_response_type_report.md),
[multivariate_response_type_report.md](multivariate_response_type_report.md),
[compositional_response_type_report.md](compositional_response_type_report.md),
[rank_choice_response_type_report.md](rank_choice_response_type_report.md),
and the interval-censored survival companion report). Unlike `nominal`, it
requires **no new `response_type` enum value at all** — it can be built as a
new `Inference*` class that operates on `response_type = "continuous"` data,
exactly the way `InferenceCountZeroInflatedNegBin` operates on
`response_type = "count"` data today without a `response_type = "zero_inflated_count"`
ever existing.

Adding **detection-limit-censored** semi-continuous support is a genuinely
different, harder project: it is a left-censoring problem, not a
point-mass-mixture problem, and the existing zero-augmented architecture
*explicitly rejects censored data* today (see `assertNoCensoring()` below).
That sub-case is **moderate to hard**, closer to the landscape survey's
original blanket rating.

## How Common Are Semi-Continuous Outcomes In Experimental Literatures?

- **Health economics / cost-effectiveness trials**: extremely common. Cost and
  utilization outcomes routinely have 20-60% of subjects at exactly $0 (no
  hospitalization, no procedure), with a right-skewed positive tail for the
  rest. The two-part model (logistic "any cost" + GLM/lognormal "cost given
  any") is close to a default analysis choice in this literature.
- **Environmental and laboratory science**: measurements below an instrument's
  detection limit are commonly recorded as zero (or as the limit itself) —
  this is the detection-limit-censored sub-case, not the structural-zero
  sub-case, and the two should not be analyzed identically even though they
  look similar in raw data.
- **Dose/utilization intensity in clinical and behavioral trials**: "amount of
  a resource used" (medication days, therapy sessions, screen time in a
  digital intervention) is naturally semi-continuous — many subjects use
  none.

Practical interpretation: this family is at least as common in applied
health-economics and utilization literatures as `count`'s own zero-inflated
variants (ZIP/ZINB) are in count literatures — and EDI already treats those
as worth first-class support. That argues for treating structural
semi-continuous outcomes the same way, not as a lower-priority edge case.

## What Already Exists

EDI already has **two proven architectural patterns** for a response that is
a mixture of a point mass plus a continuous-family component, and both are
directly relevant precedent:

### Pattern A: the count two-formula "zero-augmented" family

`InferenceCountZeroAugmentedPoissonAbstract`
([EDI/R/inference_count_zero_augmented_poisson_abstract.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_zero_augmented_poisson_abstract.R))
is the shared base class behind `InferenceCountZeroInflatedPoisson`,
`InferenceCountZeroInflatedNegBin`
([EDI/R/inference_count_zero_inflated.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_zero_inflated.R)),
and `InferenceCountHurdlePoisson`
([EDI/R/inference_count_hurdle.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_hurdle.R):16-34).
It provides, generically over the "zero-generating process" vs
"count-generating process" split:

- two independent design-matrix/formula slots — `model_formula` (conditional
  component) and `model_formula_zero` (auxiliary zero/hurdle component) —
  see `initialize()` at
  [inference_count_zero_augmented_poisson_abstract.R:30-57](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_zero_augmented_poisson_abstract.R#L30-L57)
- native Rcpp fast paths (`fast_zinb_cpp`, `fast_zero_augmented_poisson_cpp`)
  with a `glmmTMB` fallback (`fit_zero_augmented_model()`, line 718)
- a treatment-only closed-form fallback for when the full covariate-adjusted
  fit fails to converge (`fit_treatment_only_hurdle_poisson_closed_form()`,
  line 481)
- a hand-derived sandwich-variance estimator that differentiates the
  two-part log-likelihood with respect to both component parameter blocks
  (`zero_augmented_sandwich_se()`, line 403)
- full `get_likelihood_test_spec()` wiring (line 621) so LR/score/gradient
  tests and confidence-interval inversion work exactly like every other
  likelihood-backed path in the package, via `fit_null`/`score`/
  `observed_information`/`neg_loglik` hooks
- `simulate_under_lik_null()` (line 1019) for parametric-bootstrap-calibrated
  LR testing
- explicit, deliberate non-support for jackknife
  (`compute_jackknife_estimate()` at line 143 unconditionally reports
  non-estimable) because delete-one refits are too unstable for this mixture
  family — a design decision a semi-continuous class should inherit as-is
  rather than re-litigate
- an explicit `assertNoCensoring(private$any_censoring)` guard in
  `initialize()`
  ([inference_count_zero_augmented_poisson_abstract.R:46](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_zero_augmented_poisson_abstract.R#L46))
  — this is the load-bearing fact behind this report's split verdict: the
  pattern is *architected* for structural zeros, not censoring, and says so.

### Pattern B: the single-likelihood three-component mixture

`InferencePropZeroOneInflatedBetaRegr`
([EDI/R/inference_proportion_zero_one_inflated_beta.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_proportion_zero_one_inflated_beta.R):26-45)
takes a different, lighter-weight route: it inherits directly from
`InferenceAsympLikStdModCache` (the common likelihood-cache base used by most
scalar-treatment-effect MLE paths) rather than building its own two-formula
layer, and models the three-component (point mass at 0, point mass at 1,
beta interior) mixture as a single combined likelihood with an auxiliary
`model_formula_zero_one` formula (defaulting to `~ .`) for the inflation
components. This is a smaller amount of new plumbing than Pattern A, but
less reusable machinery (no shared sandwich-SE helper, no shared two-formula
matrix-building helpers) — Pattern A is the better template for a new class
specifically *because* a structural-zero-plus-continuous-tail model is
closer in shape to hurdle-Poisson (unbounded positive tail) than to
zero-one-inflated beta (bounded interior).

### Design-layer support

`Design`'s response-type validation
([EDI/R/design_abstract.R:94](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_abstract.R#L94))
already allows `continuous`. Confirmed: none of `InferenceCountZeroInflatedPoisson`,
`InferenceCountZeroInflatedNegBin`, `InferenceCountHurdlePoisson`, or
`InferenceCountHurdleNegBin` required a new `response_type` value — they all
operate on the existing `count` enum value. A structural-zero continuous
class can do the same on the existing `continuous` value. **This is the key
difference from the `nominal` report's conclusion**: there, the estimand
itself was ambiguous and required new `Design`-level storage
(`nominal_levels`) before anything else could proceed. Here, the storage is
already correct — continuous responses are already stored as plain numerics,
and a point mass at zero is just a value in that same numeric vector. No
`Design` changes are needed for the structural-zero sub-case.

## Difficulty By Layer

### 1. Design base classes: none needed (structural-zero case)

No new `response_type` value, no new storage field, no new validation. The
existing `continuous` response-type path is sufficient. This is a genuine
zero-cost layer for this feature, unlike every other family evaluated in this
batch of reports.

### 2. Most design classes: easy

Same reasoning as the nominal report: most concrete `Design` subclasses are
response-type agnostic once initialized, and since no new response type is
introduced here, there is nothing for them to adapt to.

### 3. KK21 / response-adaptive weighting designs: easy (structural-zero case)

`DesignSeqOneByOneKK21`'s response-family branches
([EDI/R/design_seq_one_by_one_KK21.R:253-281](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_seq_one_by_one_KK21.R#L253-L281))
already have a `continuous` branch. A semi-continuous outcome observed during
sequential assignment is *still a continuous response* from this layer's
point of view — no separate weighting logic is required unless response-adaptive
weighting should specifically account for the zero-inflation structure
(a refinement, not a blocker).

### 4. The new Inference class itself: moderate

This is where the real work is, but it is substantially mechanical given
Pattern A as a template:

- Choose a positive-tail continuous family for the non-zero component —
  lognormal and gamma are the standard choices in the cost/utilization
  literature (Tobit's original formulation used a censored-normal
  formulation instead, which is a materially different, harder model — see
  "The Main Conceptual Problem" below).
- Write one new C++ kernel analogous to `fast_zero_augmented_poisson_cpp`:
  a two-part negative log-likelihood (`logistic zero-process` +
  `truncated-lognormal-or-gamma positive-process`), its score, and its
  observed-information Hessian. This is new numerical work but follows an
  established shape in `EDI/src` (compare the already-implemented
  `fast_zero_augmented_poisson_cpp`/`fast_zinb_cpp` kernels), not a novel
  architecture.
- Reuse `build_component_matrix()`, `build_component_frame()`,
  `zero_augmented_sandwich_se()`'s structure (the sandwich derivation would
  need re-deriving for the new likelihood, but the surrounding
  matrix-bookkeeping code is directly reusable), and the
  `get_likelihood_test_spec()` wiring pattern verbatim.
- The `glmmTMB`-fallback path (`fit_zero_augmented_model()`) generalizes
  almost for free: `glmmTMB` supports `ziGamma(link = "log")` and a
  log-normal hurdle can be approximated via `glmmTMB`'s Gaussian family on
  `log(y)` restricted to the positive subset, mirroring
  `fit_treatment_only_hurdle_poisson_closed_form()`'s already-established
  truncated-conditional-on-positive-subset approach.

### 5. `SimulationFramework`: moderate — and not uniquely a semi-continuous gap

Confirmed via grep: **none** of the existing zero-augmented count classes
(`InferenceCountZeroInflatedPoisson`, `InferenceCountZeroInflatedNegBin`,
`InferenceCountHurdlePoisson`, `InferenceCountHurdleNegBin`) are referenced
anywhere in
[EDI/R/simulations_framework.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R) —
there is no zero-inflated/hurdle count DGP generator, no `betaT` semantics
for a two-part treatment effect, and none of the four classes appear in the
curated default inference sets. This means `SimulationFramework` support for
a new semi-continuous class is not a *semi-continuous-specific* gap — it is
a pre-existing gap shared by the entire zero-augmented-count family that
happens to also apply here. A semi-continuous DGP generator (Bernoulli zero
process + lognormal/gamma positive process, with `betaT` defined as the
treatment's effect on the conditional-positive-mean, mirroring how the new
inference class itself reports its treatment coefficient from the
conditional component only) is a moderate, well-scoped addition, but should
be scoped together with — not blocked by — adding the same support for the
existing count zero-augmented family, since the generator machinery
(Bernoulli-gated component selection) is shared.

## The Main Conceptual Problem: Point Mass vs. Censoring

For the structural-zero case, the estimand is unambiguous in a way `nominal`
and `multivariate` are not: the treatment effect is the conditional-mean
shift on the positive-tail component (exactly mirroring how the existing
zero-augmented count classes report `beta_hat_T` from the conditional
component alone — see `out$beta_hat_T = as.numeric(fit$params[2])` at
[inference_count_zero_augmented_poisson_abstract.R:838](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_count_zero_augmented_poisson_abstract.R#L838)).
No new API shape, no vector-valued output, no new summary-table
convention is needed.

The one real conceptual fork is **point mass vs. censoring**, and it changes
the model class entirely:

- **Point mass (this report's primary target)**: a zero is a *true* outcome
  value — "no cost was incurred." The generative model is a genuine mixture:
  `P(Y = 0) = p`, `P(Y = y | Y > 0) = f_positive(y)`. This is exactly what
  Pattern A already implements for counts.
- **Censoring (Tobit's original motivation, and the detection-limit case)**:
  the *true* value could be negative or a small positive number, but
  anything at or below a threshold is recorded as the threshold. The
  generative model is `Y* ~ Normal(...)`, `Y = max(Y*, 0)` (or, more
  generally, `Y = max(Y*, limit)`). This requires a censored-likelihood
  formulation (a Tobit model proper), which is architecturally a **survival-style
  censoring problem**, not a mixture problem — and the existing zero-augmented
  architecture's `assertNoCensoring()` guard confirms this package's authors
  already recognized the two do not share a code path.

Treating these as one feature would be a mistake; they should be two
separate `Inference*` classes (or families) with two separate difficulty
verdicts, which is exactly what the landscape survey's TODO-7 conflated.

## How The Existing `inference_all_*` Paths Would Respond

### `InferenceCountLikelihood` / `InferenceAsympLikStdModCache`: directly reusable

A new `InferenceContinHurdleLogNormal`-style class can inherit the same
`InferenceCountZeroAugmentedPoissonAbstract`-style abstract base (renamed or
generalized to be response-family-neutral, since its logic is not actually
count-specific except for the `za_family()`/kernel dispatch — see
"Recommended Implementation Plan" below), the same way
`InferencePropZeroOneInflatedBetaRegr` reuses `InferenceAsympLikStdModCache`
today.

### `InferenceAllSimpleMeanDiff` / `InferenceAllSimpleWilcox`: needs no new fence

Unlike `nominal`, these classes do not need a new rejection guard: they
already operate correctly (if suboptimally — ignoring the mixture structure)
on any numeric `continuous` vector, point mass at zero included. A raw mean
difference or Wilcoxon test on semi-continuous data is a valid (if
lower-powered) analysis, not a nonsensical one the way a mean-difference on
unordered category codes would be. This is a materially lower integration
risk than `nominal`.

### `InferenceSuite`: works out of the box for structural-zero case

Because no new `response_type` is introduced, `InferenceSuite`'s
compatibility-discovery mechanism
([EDI/R/inference_suite.R:110](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_suite.R#L110))
needs no changes at all for the structural-zero sub-case — a new
`continuous`-typed class is automatically discoverable the same way any new
count- or proportion-typed class already is.

## Recommended First-Wave Semi-Continuous Inference Set

1. `InferenceContinHurdleLogNormal` — logistic zero-process + log-normal
   positive-process, treatment effect on log-scale conditional mean. Closest
   analogue: `InferenceCountHurdlePoisson`.
2. `InferenceContinHurdleGamma` — logistic zero-process + Gamma(log link)
   positive-process, for right-skewed positive tails where log-normality is
   a poor fit (a common alternative in the health-economics literature).
   Closest analogue: same abstract base, different `za_family()`.
3. (Second wave, separate project) `InferenceContinTobit` — a genuine
   censored-normal Tobit model, sharing censoring-indicator plumbing with the
   survival paths rather than the zero-augmented-count plumbing. This is the
   detection-limit-censored sub-case and should not be built on the same
   abstract base as (1)/(2).

## `SimulationFramework` Implications

### Data generation

A structural-zero generator needs: draw a Bernoulli "positive" indicator per
subject (probability possibly treatment-dependent, mirroring the existing
zero-process), then draw the positive value from a log-normal or Gamma
distribution conditional on being positive. This is a direct generalization
of the existing `transform_cont_y_based_on_response_type()`
switch-based dispatch
([EDI/R/simulations_framework.R:126](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R#L126))
— conceptually smaller in scope than `nominal`'s softmax-utility generator,
since it only needs one Bernoulli gate plus the *already-existing* continuous
generator for the positive part.

### Treatment effect semantics

`betaT` is scalar and can be defined identically to how the new inference
class itself reports it: shift in the conditional-positive-mean on the log
scale, holding the zero-process's own treatment effect (if any) as a
separate, secondary nuisance parameter the DGP can also expose but need not
feature in `betaT`. No vector-valued truth is required.

### Truth and summary metrics

Because the estimand is scalar, MSE/coverage/power summaries need no
generalization beyond what already exists for every other response family.

## Recommended Implementation Plan

### Stage 1: Generalize Pattern A's abstract base (or write a lightweight sibling)

Either (a) refactor
`InferenceCountZeroAugmentedPoissonAbstract`'s truly count-specific pieces
(the `za_family()`/kernel dispatch, the Poisson/NegBin-specific score
formulas in `zero_augmented_sandwich_se()`) behind a family-pluggable
interface so a continuous sibling can share the matrix-building and
likelihood-test-spec plumbing directly, or (b) write a new, smaller abstract
base (`InferenceContinZeroAugmentedAbstract`) that copies the proven
two-formula/sandwich-SE/likelihood-test-spec shape without touching the
count-specific class. (b) is lower-risk and does not require regression
-testing the existing, working ZIP/ZINB/Hurdle-Poisson/Hurdle-NegBin paths;
(a) is more maintainable long-term. Either is a reasonable Stage 1 choice.

### Stage 2: New C++ kernel(s)

Implement `fast_zero_augmented_lognormal_cpp` (and, if pursuing both
first-wave classes, `fast_zero_augmented_gamma_cpp`) with `_score`/`_hessian`
companions, following `fast_zero_augmented_poisson_cpp`'s existing signature
shape so the R-side wiring in Stage 1 needs minimal adaptation.

### Stage 3: Concrete classes

`InferenceContinHurdleLogNormal` and (optionally in the same wave)
`InferenceContinHurdleGamma`, each a thin subclass supplying only
`za_family()`/kernel-dispatch, exactly mirroring how
`InferenceCountZeroInflatedPoisson`/`InferenceCountZeroInflatedNegBin`/
`InferenceCountHurdlePoisson` are thin subclasses of their shared abstract
base today.

### Stage 4: `SimulationFramework` support

Add the Bernoulli-gated continuous DGP generator and curated default set —
ideally scoped together with backfilling the same support for the existing
(currently DGP-less) count zero-augmented family, since the gating logic is
shared.

### Stage 5 (separate project): Detection-limit-censored Tobit

Only after Stages 1-4 are stable should a true censored-normal Tobit path be
scoped, sharing censoring-indicator plumbing with the survival module rather
than the zero-augmented-count module, and getting its own difficulty
assessment rather than inheriting this report's easier verdict.

## Bottom Line

- **Structural semi-continuous / hurdle-continuous (Tobit-like point mass at
  zero)**: **easy to moderate**. No new `response_type`, no `Design`-layer
  changes, no new fences on generic inference classes, and a directly
  reusable architectural template (`InferenceCountZeroAugmentedPoissonAbstract`)
  that has already solved every hard sub-problem (two-formula design
  matrices, sandwich SE for a two-part likelihood, likelihood-test-spec
  wiring, parametric-bootstrap null simulation, and the decision to disable
  jackknife) for a structurally identical count-family case. The genuinely
  new work is one or two C++ likelihood/score/Hessian kernels for the
  positive-tail continuous distribution, which is incremental engineering,
  not new architecture. This is the easiest of the six architectural-extension
  families evaluated across this batch of reports.
- **Detection-limit-censored (Tobit proper)**: **moderate to hard**, matching
  the landscape survey's original blanket rating. This is a censoring
  problem the existing zero-augmented pattern explicitly declines to handle
  (`assertNoCensoring()`), and should be scoped and evaluated as its own
  project sharing plumbing with the survival module instead.
- `SimulationFramework` support is moderate work, but it is not a
  semi-continuous-specific gap — the entire zero-augmented-count family
  already lacks DGP/curated-default support today, so this should be scoped
  as a shared backfill rather than a semi-continuous-only cost.
