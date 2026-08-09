# Interval-Censored Survival Response Type Report

## Scope

This report evaluates how difficult it would be to add first-class support for
**interval-censored** survival outcomes to the package, where the event time
`T` is known only to lie in an interval `(L, R]` (or, for left-censored /
current-status data, `T <= R`; for right-censored data in the ordinary sense,
`T > L`), as opposed to the package's existing `survival` response type, which
stores exactly one observed time `y` and one binary status `dead` (event
observed vs. right-censored at `y`).

This is explicitly **not** the same feature as the existing
`InferenceSurvivalDepCensTransformRegr`
([inference_survival_dep_cens_transform.R](../EDI/R/inference_survival_dep_cens_transform.R)),
which handles **dependent** right-censoring (the censoring mechanism may
correlate with the outcome) via a transformation on top of an ordinary
`(y, dead)`-shaped Cox surrogate fit
([inference_survival_dep_cens_transform.R:56-60](../EDI/R/inference_survival_dep_cens_transform.R)).
That class still consumes exactly one observed time per subject. Interval
censoring is a different structure: two time bounds per subject, with no
single "observed time" at all for subjects whose event fell strictly between
visits.

The report covers:

- what the current `survival` response storage assumes, and why it cannot
  represent an interval today
- how the two dominant survival likelihood engines already in the package
  (Cox partial likelihood, custom Weibull MLE) each interact with interval
  data
- a pragmatic, staged plan that separates a tractable parametric first step
  from a much harder semiparametric project

## Short Answer

Adding the **label** `interval_censored` (or an `interval_censoring = TRUE`
flag on the existing `survival` type) is easy.

Making a **useful subset** of interval-censored inference work — specifically
a parametric (Weibull-type) interval-censored regression path — is a
**moderate** project, because the package already has a from-scratch C++
parametric survival log-likelihood
(`fast_weibull_regression_cpp`) that can be extended with a second
log-likelihood branch, rather than needing a new numerical core built from
nothing.

Making the package's **dominant** survival inference family — Cox
proportional-hazards partial likelihood, which every one of the package's 6
non-Weibull survival paths is built on — work under interval censoring is
**hard to very hard**, because Cox's partial likelihood is fundamentally an
exact-ordering argument over a risk set, and it has no standard extension to
data where the ordering of event times is itself unknown. The standard
semiparametric fix (NPMLE / Turnbull estimator, or an EM-based
`icenReg::ic_sp`-style approach) is a materially different estimation
algorithm, not a data-format extension of `coxph.fit`.

So, mirroring `nominal_response_type_report.md`'s framing: the real
difficulty is not the enum. It is that the package's survival family is built
on two likelihood engines with very different appetites for interval data,
and the response-storage layer currently cannot even represent the input
those engines would need.

## How Common Are Interval-Censored Outcomes In Experimental Literatures?

### Clinical trials with scheduled visits

Interval censoring is the norm, not the exception, whenever an event can only
be detected at a clinic visit rather than continuously monitored:

- oncology trials measuring **progression-free survival**, where tumor
  progression is only assessed at scheduled scans — the true progression time
  lies somewhere between the last progression-free scan and the first
  progression-positive scan; see Sun (2006), *The Statistical Analysis of
  Interval-censored Failure Time Data*, Springer, the standard methodological
  reference, and its clinical motivating examples.
- HIV/AIDS cohort studies where seroconversion is only known to have occurred
  between two successive negative/positive antibody tests (the classic
  motivating example for interval-censored methods in biostatistics).
- dental and epidemiological cohort studies where a condition (e.g. caries
  onset) is only observed at periodic exams.

### Screening and diagnostic-test studies

Any study designed around periodic screening (cancer screening trials,
subclinical-disease detection) generates interval-censored onset times by
construction, since the true onset time is never directly observed.

### Reliability / engineering

Interval censoring appears wherever units are inspected on a schedule rather
than monitored continuously (e.g. periodic equipment inspection where
failure is only detected at the next inspection).

### Practical interpretation for this package

Interval censoring is one of the most common real-world departures from the
package's current `survival` assumption, because *any* experiment using
scheduled follow-up visits rather than continuous/administrative-registry
monitoring produces interval-censored, not exactly-observed, event times.
That makes it a legitimate priority gap — but the standard applied practice
in the fields above is dominated by two specific model families (parametric
AFT and semiparametric Cox-type models), which is why this report evaluates
them separately rather than proposing one generic "interval-censored"
inference path.

## What Already Exists

The package's survival family is built on exactly one response-storage shape:
a single observed time `y` and a single event/censoring indicator `dead`
(`1` = event observed at `y`, `0` = right-censored at `y`), enforced in the
`Design` base class:

- `y`, `y_original`, `dead` are initialized as parallel per-subject vectors in
  [design_abstract.R:124-128](../EDI/R/design_abstract.R).
- `add_one_subject_response()` accepts one scalar `y` and one scalar `dead`
  per subject
  ([design_abstract.R:145-194](../EDI/R/design_abstract.R)), and explicitly
  special-cases `y == 0` for `response_type == "survival"` by substituting
  `.Machine$double.eps`
  ([design_abstract.R:172-175](../EDI/R/design_abstract.R)) — a rule that
  only makes sense for a single observed time, not an interval.
- `add_all_subject_responses()` mirrors the same single-scalar-per-subject
  contract for the whole-vector case
  ([design_abstract.R:200-216](../EDI/R/design_abstract.R)), and explicitly
  rejects any censoring indicator (`dead == 0`) for non-`survival` response
  types, confirming `dead` is understood package-wide as *the* survival
  censoring flag, not a generic flag.
- The response-type enum itself is validated in exactly two places —
  [design_abstract.R:94](../EDI/R/design_abstract.R) and
  [design_abstract.R:416](../EDI/R/design_abstract.R) — both listing
  `c("continuous", "incidence", "proportion", "count", "survival", "ordinal")`.

Downstream, every survival inference path consumes `private$y` and
`private$dead` as this exact `(time, status)` pair. Two likelihood engines
dominate:

1. **Cox partial likelihood**, via `survival::coxph.fit(x = X, y =
   survival::Surv(y, dead), ...)`
   ([inference_survival_coxph.R:10-13](../EDI/R/inference_survival_coxph.R)).
   This same `Surv(y, dead)` construction, or the equivalent native
   `fast_coxph_regression_cpp(X, y, dead, ...)` kernel, is reused across
   `InferenceSurvivalCoxPH`, `InferenceSurvivalStratCoxPH`,
   `InferenceSurvivalKKStratCoxPHIVWC`
   ([inference_survival_KK_strat_cox.R:171-203](../EDI/R/inference_survival_KK_strat_cox.R)),
   `InferenceAbstractKKLWACoxIVWC`
   ([inference_survival_KK_lwa_cox_ivwc_abstract.R:200](../EDI/R/inference_survival_KK_lwa_cox_ivwc_abstract.R)),
   and `InferenceAbstractKKLWACoxOneLik`
   ([inference_survival_KK_lwa_cox_one_lik_abstract.R:138-186](../EDI/R/inference_survival_KK_lwa_cox_one_lik_abstract.R)),
   and `InferenceSurvivalDepCensTransformRegr`'s underlying surrogate fit
   ([inference_survival_dep_cens_transform.R:56-60](../EDI/R/inference_survival_dep_cens_transform.R)).
   That is **7 of 9** concrete survival-family source files.

2. **Custom parametric Weibull MLE**, via the native
   `fast_weibull_regression_cpp(y, dead, X, ...)` kernel plus matching
   `get_weibull_regression_score_cpp` / `get_weibull_regression_hessian_cpp`
   kernels wired into `get_likelihood_test_spec()`
   ([inference_survival_weibull.R:106-222](../EDI/R/inference_survival_weibull.R)).
   The same kernel is reused by `InferenceAbstractKKWeibullFrailtyIVWC`
   ([inference_survival_KK_weibull_frailty.R:184,313](../EDI/R/inference_survival_KK_weibull_frailty.R))
   and `InferenceSurvivalKKWeibullMarginal`
   ([inference_survival_KK_weibull_marginal.R:159](../EDI/R/inference_survival_KK_weibull_marginal.R)).

A third, smaller family exists outside both: `InferenceSurvivalKKClaytonCopulaIVWC`
(Clayton copula with Weibull AFT margins,
[inference_survival_KK_clayton_copula.R:4-13](../EDI/R/inference_survival_KK_clayton_copula.R))
and `InferenceAbstractKKSurvivalRankRegrIVWC` (rank regression,
[inference_survival_KK_rank_regr_ivwc_abstract.R:18-20](../EDI/R/inference_survival_KK_rank_regr_ivwc_abstract.R)).
Both are still `(y, dead)`-shaped inputs.

`SimulationFramework`'s survival data generator only ever produces ordinary
independent right-censoring: it draws a per-subject Bernoulli censoring
indicator governed by a single scalar `prob_censoring`
([simulations_framework.R:602,702](../EDI/R/simulations_framework.R)), and
the no-censoring code path explicitly sets `dead = 1` for every subject
([simulations_framework.R:3999](../EDI/R/simulations_framework.R)). There is
no concept anywhere in the generator of a second, later inspection time.

The response-adaptive `DesignSeqOneByOneKK21` weighting logic for survival
also consumes the same shape directly: `compute_weight_KK21_survival()`
builds `survival::Surv(ys_to_date, deaths_to_date)`
([design_seq_one_by_one_KK21.R:350-351](../EDI/R/design_seq_one_by_one_KK21.R))
and dispatches to the native `kk21_survival_weights_cpp(X, y, dead)` kernel
([design_seq_one_by_one_KK21.R:271-275](../EDI/R/design_seq_one_by_one_KK21.R)).

One relevant fact worth noting for later stages: R's own `survival` package
already supports interval-censored data for its **parametric** models —
`survival::survreg()` accepts `Surv(time1, time2, type = "interval2")` — but
EDI never constructs that form anywhere; every `Surv(...)` call in the
package passes exactly `(y, dead)`. So the interval-censoring capability the
package would need is not blocked by an upstream dependency gap for the
parametric case; it is blocked entirely by EDI's own response-storage and
kernel-call-site assumptions.

## Difficulty By Layer

### 1. Design base classes: moderate to hard

Unlike `nominal`, which could reuse the existing `ordinal_levels` /
factor-coercion pattern almost directly
(`nominal_response_type_report.md`, "Difficulty By Layer" §1), interval
censoring requires a genuine **schema change**, not just a new allowed enum
value plus new bookkeeping fields:

- `add_one_subject_response()` and `add_all_subject_responses()` would need a
  second time argument (e.g. `y_upper`, defaulting to `Inf` for ordinary
  right-censored subjects and to `y` itself for exactly-observed subjects, so
  existing callers with plain right-censored data continue to work
  unchanged).
- The `y == 0` → `.Machine$double.eps` substitution at
  [design_abstract.R:172-175](../EDI/R/design_abstract.R) needs an
  interval-aware analogue (what does a zero-width interval starting at 0
  mean? Left-censored/current-status data, most likely — that is its own
  special case worth naming explicitly rather than silently coercing).
- `assert_y()` needs an interval-consistency check (`y_lower <= y_upper`,
  and a coherent policy for `Inf` upper bounds representing ordinary
  right-censoring, `0` lower bounds representing left-censoring, and
  `y_lower == y_upper` representing an exactly observed event — this is the
  standard 4-type interval-censoring taxonomy used by `survival::Surv(...,
  type = "interval2")` and should be adopted rather than reinvented).
- Every downstream site that currently reads `private$y` and `private$dead`
  as *the* survival representation (all 9 files above) would need to either
  read a new pair of fields, or the package would need to keep `(y, dead)`
  as an ordinary-right-censoring-only representation and introduce a
  **parallel** `(y_lower, y_upper)` pair used only by interval-aware paths —
  the latter is safer (zero risk of silently breaking the 7+ existing
  Cox/Weibull call sites) and is the approach this report recommends.

### 2. Most design classes: easy

As with `nominal`, design/randomization classes that only pass responses
through (assignments, covariates, sequential arrival) do not interpret `y`
or `dead` themselves — they would be unaffected by adding an optional second
time field, provided the field defaults sensibly for the ~100% of use cases
that remain ordinary right-censored data.

### 3. KK21 / response-adaptive weighting designs: hard

`compute_weight_KK21_survival()` builds a `survival::Surv(y, dead)` object
and a matching native C++ weighting kernel
([design_seq_one_by_one_KK21.R:350-351,271-275](../EDI/R/design_seq_one_by_one_KK21.R)).
Response-adaptive weighting under interval censoring would need either a
genuinely interval-aware weighting scheme (nontrivial — even the applied
interval-censoring literature does not have a single standard "on-the-fly"
covariate-weighting recipe the way it does for right-censored survival), or
an explicit documented approximation (e.g. weight using the interval
midpoint or upper bound as a plug-in "observed time," which is a known but
biased shortcut). This is a harder design-side item than for `nominal`,
because there is no clean fallback that avoids inventing new statistical
methodology.

### 4. `SimulationFramework`: hard, but narrower than for `nominal` or `multivariate`

Unlike `nominal`'s truth-object problem, interval censoring does **not**
change the treatment-effect estimand — `betaT` remains whatever scalar
survival effect (log-hazard-ratio or log-time-ratio) the package already
uses. What is genuinely new work:

- `prob_censoring`
  ([simulations_framework.R:602,702](../EDI/R/simulations_framework.R)) is a
  single independent-right-censoring-probability generator; an
  interval-censoring generator instead needs a **visit-schedule** parameter
  (e.g. inspection times at a fixed or random cadence) and must report,
  for each subject, the last pre-event visit and first post-event visit
  rather than an exact time.
- The no-censoring shortcut that hard-codes `dead = 1` for every subject
  ([simulations_framework.R:3999](../EDI/R/simulations_framework.R)) and the
  mean-diff exclusion rule keyed on `prob_censoring > 0`
  ([simulations_framework.R:4137](../EDI/R/simulations_framework.R)) would
  both need an interval-censoring-aware branch.
- Truth/coverage/MSE summary logic does not need new *semantics* (the truth
  is still a scalar `betaT`), only a new *data-generating path* feeding into
  the same downstream summary code — this is meaningfully less architectural
  work than the `nominal` or `multivariate` cases, where the estimand itself
  is undefined until a design choice is made.

## The Main Conceptual Problem: Cox Partial Likelihood Does Not Extend To Interval Data

For every other response family in the package, the "difficulty" question is
mostly about plumbing. For interval-censored survival, there is a genuine
**statistical** blocker sitting underneath the plumbing question, and it
falls on the package's single largest survival sub-family (Cox partial
likelihood, 7 of 9 survival source files per "What Already Exists" above):

- Cox's partial likelihood is constructed from the *ordering* of observed
  event times within risk sets. It requires knowing, at each observed event
  time, exactly which other subjects were still at risk and had not yet
  failed. Interval-censored data does not provide an exact ordering of event
  times — two subjects' true event times may be unorderable given only their
  overlapping intervals.
- There is no drop-in interval-censored analogue of `coxph.fit()` in the
  `survival` package the way there is for the parametric AFT case
  (`survreg()` already accepts `Surv(..., type = "interval2")`, as noted
  above). The standard semiparametric approaches — Turnbull's NPMLE, or an
  EM-based proportional-hazards fit as implemented in the CRAN package
  `icenReg` (`ic_sp()`) — are different estimation algorithms with different
  convergence behavior, not a data-format change to the existing
  `coxph.fit()` call.
- This means the 7 Cox-based survival paths cannot be extended to interval
  data by touching the response-storage layer alone; each would need either
  a new semiparametric estimation engine (large, genuinely new statistical
  software, likely its own C++ kernel and its own optimizer convergence
  story) or a documented approximation (e.g. treating interval midpoints as
  exact times and accepting the resulting bias — a common but methodologically
  weak shortcut that should be labeled clearly as approximate, not silently
  offered as exact inference).

By contrast, the **parametric** Weibull path does not have this problem:
maximum likelihood for an interval-censored observation is simply
`log(S(L; theta) - S(R; theta))` in place of the usual
`log(f(y; theta))` / `log(S(y; theta))` terms, which is a closed-form,
differentiable modification to an existing from-scratch C++ log-likelihood
(`fast_weibull_regression_cpp` plus its paired score/Hessian kernels,
[inference_survival_weibull.R:106-222](../EDI/R/inference_survival_weibull.R)).
This is why this report gives sharply different verdicts for the parametric
vs. semiparametric cases below, rather than one blanket "hard" rating for
"interval-censored survival" as a whole.

## How The Existing `inference_all_*` Paths Would Respond

### `InferenceAsympLikStdModCache` (Weibull's parent): reusable, with a new likelihood branch

`InferenceSurvivalWeibull` and the frailty/marginal Weibull variants already
inherit the standard likelihood-caching abstraction used across the package
(the same one `nominal_response_type_report.md` identifies as reusable for
scalar likelihood-based models). Extending the Weibull kernel's
log-likelihood with an interval-censoring branch requires no change to this
abstraction — it already expects a scalar coefficient vector, a
`get_likelihood_test_spec()` implementation, and score/Hessian callbacks,
all of which the interval-censored Weibull path would still provide in the
same shape.

### Cox-based paths (`InferenceSurvivalCoxPH` and its 6 relatives): cannot accept interval data without a new estimation engine

None of these paths can be made interval-aware by changing what they are fed
— they are constructed directly around `coxph.fit()` / partial likelihood.
Making them "respond" to interval data correctly means either (a) explicitly
rejecting interval-censored designs (the safe, cheap option, matching how
`nominal_response_type_report.md` recommends generic paths explicitly reject
response types they cannot honor), or (b) building a wholly new
semiparametric estimation path that happens to share the same class-naming
convention.

### `InferenceSurvivalDepCensTransformRegr`: should reject, and should say why

Because this class's surrogate fit is itself Cox-based
([inference_survival_dep_cens_transform.R:56-60](../EDI/R/inference_survival_dep_cens_transform.R)),
it inherits the same Cox limitation. Its docstring/roxygen should gain an
explicit note that it addresses *dependent censoring*, not *interval
censoring*, since the two are easy to conflate from the name alone — this
was exactly the ambiguity flagged in TODO-8 of
`response_types_landscape_report.md`.

### Randomization / bootstrap base classes: mostly neutral, same as for `nominal`

The generic bootstrap/randomization machinery does not care about the
internal shape of the likelihood being resampled; it only needs a scalar
statistic and a way to refit under resampled weights. So these base classes
are not a blocker for either the parametric or semiparametric interval
paths, mirroring the finding in `nominal_response_type_report.md`'s
equivalent section.

### `InferenceSuite`: needs explicit interval-awareness fences

Exactly as with `nominal`, `InferenceSuite` discovers applicability by
instantiation and error-catching
([inference_suite.R:110](../EDI/R/inference_suite.R)). If a
`survival`-typed design carries an interval-censored response and a
Cox-based class is not taught to reject it explicitly, `InferenceSuite` may
silently include a Cox path that will either error unpredictably deep inside
`coxph.fit()` or (worse) silently coerce the interval midpoint into a single
`y` and produce a plausible-looking but methodologically wrong answer. This
argues strongly for making interval-censored responses a **distinguishable,
explicit** state on the design object (not just an implicit consequence of
`y_upper` being set) so every existing survival class can assert against it
up front.

## Package Areas That Need Explicit Interval-Censoring Fences

Before implementing any interval-aware inference path, the following should
explicitly reject (or, at minimum, loudly warn on) an interval-censored
design:

- `InferenceSurvivalCoxPH`, `InferenceSurvivalStratCoxPH`,
  `InferenceSurvivalKKStratCoxPHIVWC`,
  `InferenceAbstractKKLWACoxIVWC` / `InferenceAbstractKKLWACoxOneLik`
  (all Cox-based)
- `InferenceSurvivalDepCensTransformRegr` (Cox-surrogate-based; also needs
  the docstring clarification above)
- `InferenceSurvivalKMDiff` and `InferenceSurvivalRMST` — both are
  Kaplan-Meier-based, and the KM estimator has its own interval-censored
  analogue (Turnbull), so these should reject rather than silently treat
  interval bounds as exact times
- `InferenceAbstractKKSurvivalRankRegrIVWC` (rank regression assumes an
  orderable exact-or-right-censored time, the same problem as Cox)
- `InferenceSurvivalKKClaytonCopulaIVWC` (Weibull-margin copula — plausibly
  extensible later using the same interval-likelihood modification as plain
  Weibull, but should reject until that extension is actually built)

## Common New Inference Paths To Implement

### 1. Interval-censored parametric Weibull regression

Extend `fast_weibull_regression_cpp` (and its score/Hessian kernels) with an
interval-censoring likelihood branch: for interval-censored subjects,
contribute `log(S(L; theta) - S(R; theta))`; for exactly-observed subjects
(`y_lower == y_upper`), contribute the ordinary density term unchanged; for
ordinary right-censored subjects (`y_upper == Inf`), contribute the ordinary
survival-function term unchanged. This is a strict generalization of the
existing kernel, so ordinary right-censored callers are unaffected.

Difficulty: **moderate**

Why: the numerical core, optimizer, and `get_likelihood_test_spec()`
plumbing already exist; only the per-subject log-likelihood term needs a new
case, and the Weibull survival/density functions needed for the interval
term are already implemented for the ordinary case.

### 2. Interval-censored parametric Weibull frailty / KK-marginal variants

Once the base kernel above exists, `InferenceAbstractKKWeibullFrailtyIVWC`
and `InferenceSurvivalKKWeibullMarginal`, which already call the same
`fast_weibull_regression_cpp` kernel, would need comparatively little
additional work to accept the new likelihood branch.

Difficulty: **easy, once path #1 exists**

### 3. Semiparametric interval-censored proportional hazards (Cox-analogue)

Implement a genuinely new estimation engine mirroring `icenReg::ic_sp()` /
Turnbull-style NPMLE-based proportional hazards, as a new class (e.g.
`InferenceSurvivalIntervalCensCoxNPMLE`) rather than modifying the existing
`InferenceSurvivalCoxPH`.

Difficulty: **hard**

Why: this is a new estimation algorithm (EM over a discretized baseline
survival function, or a spline-based semiparametric MLE), with its own
convergence and standard-error story, not an extension of `coxph.fit()`.
This is comparable in scope to the package's existing GEE or KK-GLMM
projects — a new statistical method, not a data-shape accommodation.

### 4. Current-status (case-I interval-censoring) special case

Where every subject is inspected exactly once (`y_lower = 0` or
`y_lower` = last known event-free time, `y_upper` = the single inspection
time), a simpler isotonic-regression-based nonparametric estimator (the
"pool-adjacent-violators"-based current-status MLE) exists and is
substantially simpler to implement than general interval censoring.

Difficulty: **moderate**

Why: this is a well-known, narrower special case with simpler estimation
theory; worth calling out separately because it may be the more common case
in some experimental designs (e.g. single end-of-study assessment) and does
not require the general interval machinery.

## Recommended First-Wave Set

Mirroring the staged approach in `nominal_response_type_report.md`, the best
first wave avoids committing to the hardest (semiparametric) case up front:

1. Admit interval-censored responses safely at the `Design` layer via a new
   optional `(y_lower, y_upper)` pair, defaulting to ordinary
   `(y, dead)`-equivalent behavior for all existing callers.
2. Add explicit rejection guards to every Cox-based and KM-based survival
   class (see "Package Areas That Need Explicit Interval-Censoring Fences").
3. Implement the interval-censored parametric Weibull path (New Path #1
   above) as the first genuinely working interval-censored inference method.
4. Only after the parametric case is stable, evaluate whether to commission
   the semiparametric Cox-analogue (New Path #3) as its own project.

Concretely:

- `InferenceSurvivalWeibullIntervalCens` (or an `interval_censoring`
  constructor flag on the existing `InferenceSurvivalWeibull`)
- explicit `stop()`/incompatibility guards added to the 7 Cox/KM-based
  classes listed above

## `SimulationFramework` Implications

### Data generation

`prob_censoring`'s independent-Bernoulli right-censoring generator
([simulations_framework.R:602,702](../EDI/R/simulations_framework.R)) needs
a sibling generator that instead draws a visit schedule (e.g. fixed-cadence
inspection times, or a Poisson-process arrival of inspection visits) and
reports, per subject, the bracketing interval around the latent event time
rather than the exact time.

### Treatment effect semantics

Unlike `nominal` or `multivariate`, **no new estimand semantics are needed**
— `betaT` remains the same scalar log-hazard/log-time-ratio effect already
used for `survival`. This is the main reason interval censoring is a
narrower `SimulationFramework` project than either of those two families.

### Truth and summary metrics

Unaffected — MSE/coverage/power summaries operate on the same scalar `betaT`
as today; only the *data-generating* step changes.

## Recommended Implementation Plan

### Stage 1: Admit the shape safely

Add an optional `(y_lower, y_upper)` pair to `Design`'s response storage,
with `assert_y()` enforcing the 4-type interval-censoring taxonomy
(exact/right/left/interval) used by `survival::Surv(..., type =
"interval2")`. Leave `(y, dead)` untouched for all existing right-censored
callers.

### Stage 2: Fence off incompatible paths

Add explicit rejection guards to the 7 Cox-based and 2 KM-based survival
classes identified above, with an error message naming the specific
incompatibility (partial likelihood requires exact event ordering).

### Stage 3: Implement the parametric interval-censored Weibull path

Extend `fast_weibull_regression_cpp` with the interval-likelihood branch;
wire it into a new class or constructor flag; extend
`get_weibull_regression_score_cpp` / `_hessian_cpp` to match.

### Stage 4: Extend `SimulationFramework` for interval data generation

Add the visit-schedule generator; no new truth/estimand logic is required.

### Stage 5 (separate, larger project): semiparametric interval-censored Cox analogue

Only after Stage 3 is stable and adopted, evaluate commissioning a dedicated
project for an NPMLE/EM-based proportional-hazards fit under interval
censoring, comparable in scope to the package's other from-scratch
estimation-engine projects.

## Bottom Line

Interval-censored survival is **not one difficulty level** — it splits
cleanly along the package's own existing survival-engine boundary:

- **Parametric (Weibull) interval-censored regression: moderate.** The
  numerical core, optimizer plumbing, and likelihood-test-spec machinery
  already exist; only the per-subject log-likelihood term needs a new
  branch.
- **Current-status (case-I) special case: moderate**, and worth
  implementing as its own simpler path rather than folding into the general
  interval machinery.
- **Semiparametric Cox-analogue interval-censored regression: hard.** Cox
  partial likelihood has no data-shape extension to interval data; this
  requires a genuinely new estimation algorithm (NPMLE/EM), affecting 7 of
  the package's 9 survival source files if full parity with the existing
  Cox family is the goal.
- **`Design`-layer response storage: moderate to hard**, because — unlike
  `nominal`, which could mostly reuse the existing ordinal factor-storage
  pattern — interval censoring needs an actual second-time-bound field, a
  real schema change rather than additive bookkeeping.
- **`SimulationFramework`: hard for the data-generation mechanism, but
  narrower than for `nominal`/`multivariate`**, because the treatment-effect
  estimand itself does not change — only the observation mechanism does.

So the pragmatic recommendation, consistent with TODO-8's framing in
`response_types_landscape_report.md`, is to treat this as a **later
survival-subfamily project**, staged so that the tractable parametric case
ships first and the genuinely hard semiparametric case is scoped and
commissioned separately, rather than as a single monolithic "add interval
censoring" effort.
