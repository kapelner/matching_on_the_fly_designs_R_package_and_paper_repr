# Longitudinal / Repeated-Measures Response Type Report

## Scope

This report evaluates how difficult it would be to add first-class support
for **longitudinal / repeated-measures outcomes** to the package — outcomes
where each subject contributes **more than one** observation over time
(panel-survey waves, repeated clinic-visit biomarker draws, weekly
sentiment measurements, etc.), rather than the single scalar (or
single-censored-time) outcome every current `response_type` assumes.

This is not a request from
[response_types_landscape_report.md](references/response_types_landscape_report.md)'s
"Additional Response Families Worth Noting" section (TODO-3 in that
document's Implementation TODOs), motivated by the "Political
experiments"/"Sociology" panel-survey field sections and the "Clinical
trials" repeated-biomarker field section. That document rated this `hard`
qualitatively; this report supplies the detailed architectural assessment.

The report covers:

- what "repeated measures" actually requires structurally, distinct from
  every existing response type
- whether the package's existing clustered-correlation machinery (KK GEE,
  KK GLMM) is reusable scaffolding — this turns out to be the central
  finding
- difficulty by architectural layer
- a pragmatic staged implementation plan

## Short Answer

Adding repeated-measures support is **hard**, but for a specific and
narrower reason than "nominal" was hard. Nominal was hard because the
package had no natural *scalar estimand* for an unordered category. Here,
the estimand is usually fine (a treatment coefficient in a marginal or
mixed model) — what's hard is that **every current response type assumes
exactly one observation per subject**, and that assumption is not a
validation nicety sitting at the edge of the code. It is load-bearing in
the `Design` base class's core response container itself.

The good news, and the main finding of this report: the package **already
has real within-cluster-correlation machinery** — GEE with a working
correlation structure (`corstr = "exchangeable"`,
[inference_mixin_kk_gee_shared.R:426](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_mixin_kk_gee_shared.R:426))
and GLMM with a subject-level random intercept
(`(1 | group_id)`,
[inference_continuous_KK_glmm.R:5](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_continuous_KK_glmm.R:5)).
Both exist today to handle correlation **induced by the KK matched-pair
design** (correlation *across different subjects* who share a match
group), not correlation from **repeated observations on the same
subject**. But the numerical engines — GEE with an exchangeable or AR(1)
working correlation, a mixed model with a random intercept per cluster —
are exactly the tools repeated-measures analysis needs. The clustering
*unit* would change (subject-over-time instead of matched-pair), but the
fitting machinery would not need to be reinvented.

So the real cost is almost entirely in the **response-storage layer**
(`Design$y`, one value per subject index `t`) and in
`SimulationFramework`'s single-scalar-per-subject data-generating process
— not in writing a new estimator from scratch.

## How Common Are Repeated-Measures Outcomes In Experimental Literatures?

### Clinical trials

Repeated biomarker draws, longitudinal quality-of-life instruments, and
multi-visit symptom scores are a mainstream primary- or secondary-endpoint
family in trials with more than one follow-up visit — arguably more common
than single-timepoint continuous endpoints once a trial has any scheduled
follow-up structure beyond a single primary endpoint visit. Mixed-effects
models for repeated measures (MMRM) are the de facto standard analysis
approach for continuous longitudinal endpoints in confirmatory trials.

### Political science / sociology field experiments

Panel surveys — the same respondent interviewed at multiple waves relative
to a treatment (e.g. a mailer, an ad exposure, an information
intervention) — are a standard design in field experiments measuring
attitude or behavior change over time, exactly the "Political
experiments"/"Sociology" motivation cited in
[response_types_landscape_report.md](references/response_types_landscape_report.md).
The treatment effect of interest is often the *trajectory* (slope) or the
*persistence* of an effect across waves, not just a single post-treatment
mean difference — a genuinely different estimand shape than anything the
package currently computes.

### Practical interpretation for this package

Unlike nominal outcomes (present but not dominant), repeated-measures
designs are extremely common **whenever a design has more than one
scheduled follow-up**, which is itself a design-level feature EDI does not
currently model at all — every `Design` subclass assumes exactly one
outcome collection event per subject
(`add_one_subject_response(t, y, dead = 1)`,
[design_abstract.R:145](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_abstract.R:145)).
So this feature is not "add a response type" so much as "add a second
design dimension" (time), which is a more structural change than any of
the six existing response types individually represent.

## What Already Exists

Two pieces of existing machinery are directly relevant, both built to
handle **within-cluster correlation**, just for a different clustering
unit (KK matched-pair, not subject-over-time):

### GEE with a working correlation structure

`InferenceMixinKKGEEShared`
([inference_mixin_kk_gee_shared.R:18](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_mixin_kk_gee_shared.R:18))
is spliced into every KK-GEE daughter class (`InferenceCountKKGEE`,
`InferenceOrdinalKKCLMMAbstract`'s GEE cousins, etc.) and fits with
`corstr = "exchangeable"` at
[inference_mixin_kk_gee_shared.R:426](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_mixin_kk_gee_shared.R:426)
and
[:476](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_mixin_kk_gee_shared.R:476).
GEE with an exchangeable, AR(1), or unstructured working correlation is
the textbook marginal-model approach to repeated measures. The cluster
variable it currently groups on is the KK match id, not a subject-time
panel id, but the fitting and sandwich-variance code does not care what
the cluster variable *means* — it only needs a cluster index vector and a
long-format data frame. That is architecturally the single biggest
head start this project would have.

### GLMM with a subject-level random intercept

`InferenceContinKKGLMM`
([inference_continuous_KK_glmm.R:1](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_continuous_KK_glmm.R:1))
fits `(1 | group_id)` where `group_id` is literally the KK match-id vector
(`group_id = m_vec`,
[inference_continuous_KK_glmm.R:160](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_continuous_KK_glmm.R:160)).
A repeated-measures random-intercept (or random-intercept-and-slope) model
is structurally the same fitting problem with `group_id = subject_id`
instead. The Rcpp/L-BFGS Gaussian LMM kernel underneath this class is very
likely reusable near-verbatim; what changes is only how `group_id` gets
constructed upstream, which is a `Design`/response-storage problem, not a
fitting-kernel problem.

Neither of these is a public, response-type-agnostic "give me a clustering
variable and a working-correlation family" utility today — both are
private to their KK-specific classes — so extracting a shared,
cluster-index-driven GEE/GLMM core is itself part of the work, but it is
refactoring existing correct code, not inventing new numerics.

## Difficulty By Layer

### 1. `Design` base class core response storage: hard

This is the layer that makes repeated measures fundamentally different
from every prior response-type addition. `private$y` is a flat vector
indexed by subject order `t`
([design_abstract.R:145](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_abstract.R:145),
`add_one_subject_response = function(t, y, dead = 1)`), and
`add_all_subject_responses(ys, deads = NULL)`
([design_abstract.R:200](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_abstract.R:200))
takes one scalar per subject. Every downstream consumer of `Design`
(inference classes, `SimulationFramework`, summary tables) reads `y` as
"the outcome," singular, per subject.

Supporting repeated measures means one of:

- extending `Design` to store `y` as a **ragged list-of-vectors** (or a
  long-format `(subject_id, time, y)` triple store) alongside the existing
  scalar path, with `response_type`-conditional dispatch — every method
  that currently does `private$y[t]` needs a repeated-measures-aware
  counterpart, and
- deciding whether the **design/randomization** side even changes (KK
  designs adapt covariate weights based on the *observed response*; for a
  longitudinal outcome, "the observed response so far" is ambiguous mid-panel
  — is it the latest wave? the running mean? This is a genuine new design
  decision, not just a storage change)

This is real, non-mechanical architectural work, and it is the layer that
makes the overall verdict `hard` rather than `moderate`.

### 2. Most design classes: easy once `Design` supports it

As with nominal, most concrete `Design` subclasses (fixed randomization,
blocked randomization, simple sequential designs) do not model the
response directly — they consume assignments and covariates. Once the base
class supports a ragged/long-format response, most of these subclasses
would need only pass-through changes.

### 3. KK21 / response-adaptive weighting designs: hard

`DesignSeqOneByOneKK21` branches explicitly on `response_type` to compute
adaptive covariate weights from the observed response
([design_seq_one_by_one_KK21.R:259](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_seq_one_by_one_KK21.R:259)
through
[:292](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_seq_one_by_one_KK21.R:292)).
This logic assumes one committed response value is available to weight on
by the time the next subject is assigned. A subject who is still mid-panel
(has some but not all waves observed) breaks that assumption outright —
this is a harder version of the same problem nominal caused for KK21
weighting, because there isn't even a well-defined "the response" to
compute a weight from until the panel is complete or a specific wave/
summary is designated. **This should probably be explicitly out of scope
for a first version**: repeated-measures outcomes could be supported only
on non-response-adaptive designs initially, with KK21-style adaptive
weighting on panel data deferred indefinitely.

### 4. `SimulationFramework`: hard

`SimulationFramework` generates one `response_type`-conditional scalar per
subject via
`transform_cont_y_based_on_response_type()`
([simulations_framework.R:126](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R:126),
`switch(response_type, ...)` at
[:136](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R:136))
and defines the scalar "true effect" used for MSE/coverage/power per
response type in a second `switch(private$current_response_type, ...)` at
[simulations_framework.R:3913](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R:3913).
Both would need a genuinely new branch: generating `T` correlated draws
per subject (e.g. an AR(1) or compound-symmetric latent process across
waves) rather than one draw, and defining what the "true treatment effect"
scalar means when the DGP has both a between-subject treatment effect and
a within-subject time trend. This mirrors the same category of cost the
nominal report found for `SimulationFramework`, and for the same root
reason: the framework's data-generation and truth-definition code is
built around exactly one scalar outcome per subject.

## The Main Conceptual Problem: What Is The Repeated-Measures Estimand?

Unlike nominal (where the problem was *no* natural scalar estimand),
repeated measures has several well-established scalar estimands to choose
from, but the package must pick one (or support more than one) explicitly:

- **marginal treatment effect on the mean outcome, averaged over time**
  (a GEE population-averaged effect) — closest in spirit to what every
  existing response type already reports
- **treatment-by-time interaction** (does the treatment effect grow,
  shrink, or stay flat across waves?) — the estimand panel-survey and MMRM
  literature usually actually cares about, and one the package has no
  precedent for (every current estimand is a single scalar with no time
  dimension)
- **subject-specific (conditional) treatment effect** from a mixed model
  — different from the marginal GEE effect whenever the link function is
  non-identity (logistic, log-linear repeated measures), same
  marginal-vs-conditional distinction that already exists between GEE and
  GLMM elsewhere in the package

A pragmatic first version should target the marginal, averaged-over-time
effect specifically because it is a single scalar with the same shape as
every existing `beta_hat_T`, and because the GEE machinery to compute it
already exists in a directly adaptable form (§ What Already Exists).
Treatment-by-time interaction is a legitimate and valuable second-wave
target but is vector-valued (one coefficient per wave or per polynomial
time term) and would need the same kind of "does the package support
vector-valued estimands" decision the nominal report flagged for
multinomial logit — this project should not try to solve that generically
on day one either.

## How The Existing `inference_all_*` Paths Would Respond

### Generic scalar-mean/rank paths (`InferenceAllSimpleMeanDiff`, `InferenceAllSimpleWilcox`): should reject or need a pre-aggregation step

These classes have no concept of repeated observations per subject; they
would need to either explicitly reject a repeated-measures response type,
or (more usefully) accept it only after the caller supplies a
per-subject-summary reduction (e.g. "test on each subject's final-wave
value" or "test on each subject's within-subject mean") — which is really
a *different*, simpler feature (subject-level summarization) layered in
front of the existing scalar machinery, not new machinery on the
mean-diff/Wilcoxon side itself.

### `InferenceAsympLikStdModCache`: reusable only for a marginal (GEE-shaped) scalar effect

This abstract class assumes `generate_mod()` returns a fitted object with
`b[2]` as a scalar treatment effect
([inference_all_abstract_asymp_lik_std_mod_cache.R:6](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_asymp_lik_std_mod_cache.R:6)).
A marginal GEE-style repeated-measures fit produces exactly that shape (a
scalar population-averaged treatment coefficient), so this layer is
reusable for the recommended first estimand above — it is not reusable
for a treatment-by-time vector effect without the same generalization the
nominal report already flagged as a shared, deferred problem.

### `InferenceMixinKKGEEShared` / `InferenceContinKKGLMM`: the direct reuse targets

As detailed above, these are the two classes worth generalizing rather
than writing repeated-measures fitting code from scratch. The concrete
engineering task is extracting their GEE-fit-with-working-correlation and
GLMM-fit-with-random-intercept logic into a response-type-agnostic mixin
parameterized by an arbitrary cluster-index vector, then supplying
`subject_id` (instead of KK match id) as that cluster index for a new
repeated-measures inference class.

### `InferenceSuite`: adapts once repeated-measures classes are response-type-gated

Same conclusion as the nominal report: `InferenceSuite` discovers
compatibility by construction-time rejection
([inference_suite.R:110](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_suite.R:110)),
so this is not a blocker as long as new repeated-measures classes assert
their response type explicitly and existing scalar classes reject it.

## Package Areas That Need Explicit Repeated-Measures Fences

Once `Design` can store a ragged/long-format response, every existing
inference class that reads `private$y` as a flat per-subject scalar vector
needs to either explicitly reject the new response type or explicitly
opt into a defined pre-aggregation (e.g. "use the final wave"). Given the
scale of the existing class roster (dozens of concrete `Inference*`
classes across `continuous`/`incidence`/`count`/`proportion`/`survival`/
`ordinal`), auditing and gating all of them is itself a large mechanical
task, structurally identical to the fencing work the nominal report
already identified as necessary there — this cost recurs for every new
response family that changes the shape of `y`, and is not unique to
repeated measures.

## Recommended First-Wave Repeated-Measures Inference Set

1. **`InferenceLongitudinalGEE`** — marginal (population-averaged)
   treatment effect via GEE with a caller-selectable working correlation
   (`exchangeable`, `ar1`, `independence`), generalizing
   `InferenceMixinKKGEEShared`'s existing fit code
   ([inference_mixin_kk_gee_shared.R:426](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_mixin_kk_gee_shared.R:426))
   to an arbitrary subject-time cluster index. This should be the
   flagship first path — it reuses the most existing machinery and
   produces the estimand shape (`beta_hat_T` scalar) every downstream
   package feature (summary tables, `InferenceAsympLikStdModCache`,
   `InferenceSuite`) already expects.
2. **`InferenceLongitudinalGLMM`** — subject-specific treatment effect via
   a random-intercept mixed model, generalizing `InferenceContinKKGLMM`'s
   `(1 | group_id)` fit
   ([inference_continuous_KK_glmm.R:160](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_continuous_KK_glmm.R:160))
   to `group_id = subject_id`.
3. Defer treatment-by-time interaction, random slopes, and any
   response-adaptive-design integration (KK21 weighting on panel data) to
   a second wave, for the same reason the nominal report deferred
   full multinomial support: they require a vector-valued or
   design-time-sequencing decision the package does not currently have an
   answer for.

## `SimulationFramework` Implications

### Data generation

A new branch is needed that generates `T` within-subject-correlated draws
(e.g. compound-symmetric or AR(1) latent noise across waves) instead of
one draw per subject, extending
`transform_cont_y_based_on_response_type()`
([simulations_framework.R:126](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R:126)).
This is not a one-line `switch` branch addition — it changes the shape of
what the generator returns (a `subject x wave` matrix, not a vector).

### Treatment effect semantics and truth

`betaT` and the true-effect computation at
[simulations_framework.R:3913](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R:3913)
are both currently scalar-per-response-type. For the recommended
marginal-GEE first estimand, `betaT` can stay scalar (a
population-averaged shift applied at every wave), which keeps this
tractable; a treatment-by-time-varying `betaT` (the more scientifically
interesting case) would need `SimulationFramework`'s parameter grid and
truth machinery generalized to vector-valued `betaT`, which is out of
scope for a first version.

## Recommended Implementation Plan

### Stage 1: Extract a response-type-agnostic clustered-fit core

Refactor the GEE-fit-with-working-correlation logic in
`InferenceMixinKKGEEShared`
([inference_mixin_kk_gee_shared.R:18](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_mixin_kk_gee_shared.R:18))
and the GLMM-fit-with-random-intercept logic in `InferenceContinKKGLMM`
([inference_continuous_KK_glmm.R:1](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_continuous_KK_glmm.R:1))
into shared mixins parameterized by an arbitrary cluster-index vector,
without changing behavior for the existing KK match-id use case. This is
low-risk, mechanical refactoring that pays for itself even before
repeated measures exists, because it is currently duplicated,
KK-specific-named logic.

### Stage 2: Add a ragged/long-format response container to `Design`

Add `response_type = "longitudinal"` (or a `panel = TRUE` modifier on an
existing type) with a `(subject_id, wave, y)`-shaped storage path,
explicit fencing on `DesignSeqOneByOneKK21`'s response-adaptive weighting
(reject or require a specified summary-wave), and explicit rejection
guards on every existing scalar-`y` inference class.

### Stage 3: Add `InferenceLongitudinalGEE`

Wire Stage 1's generalized GEE mixin to Stage 2's storage, targeting the
marginal population-averaged treatment effect as the estimand.

### Stage 4: Add `InferenceLongitudinalGLMM`

Wire Stage 1's generalized GLMM mixin similarly, targeting the
subject-specific (conditional) treatment effect.

### Stage 5: Extend `SimulationFramework`

Only after Stages 2-4 land should the DGP and truth-definition layers be
extended, mirroring the nominal report's own sequencing rationale:
building the simulation harness before the estimand and storage model are
settled would mean rebuilding it once they change.

## Bottom Line

Adding longitudinal / repeated-measures support is **hard**, for a
different reason than nominal was hard: the blocker here is not "no
scalar estimand" (a marginal GEE effect is a perfectly good,
already-precedented scalar estimand) — it is that **every existing
response type, and the `Design` base class's core response container
itself, hard-assumes exactly one observation per subject**. That
assumption is structural, not a validation nicety, and touching it has
ripple effects into KK21 adaptive weighting and `SimulationFramework`'s
data generation, both of which are themselves rated `hard` in isolation.

The mitigating finding is that the **numerical fitting engines already
exist** in a directly adaptable form: `InferenceMixinKKGEEShared`'s
exchangeable-working-correlation GEE fit and `InferenceContinKKGLMM`'s
random-intercept GLMM fit were both built to handle KK-matched-pair
clustering, and generalizing their cluster-index parameter to a
subject-over-time panel id is a comparatively modest, low-risk refactor.
So the pragmatic recommendation is:

1. extract the existing GEE/GLMM clustered-fit logic into
   response-type-agnostic, cluster-index-parameterized mixins first (pays
   for itself immediately as a cleanup, independent of this feature)
2. add a ragged/long-format response container to `Design`, fencing off
   response-adaptive (KK21) designs and every existing scalar-`y`
   inference class explicitly
3. ship a marginal-GEE repeated-measures inference class first, since it
   reuses the most existing machinery and produces a scalar estimand every
   downstream package feature already expects
4. add a random-intercept GLMM repeated-measures class second
5. postpone treatment-by-time interaction, random slopes, and
   response-adaptive-design integration until the vector-valued-estimand
   question (already deferred by the nominal report for the same
   underlying reason) is answered generally
6. postpone `SimulationFramework` support until the storage model and
   estimand are settled
