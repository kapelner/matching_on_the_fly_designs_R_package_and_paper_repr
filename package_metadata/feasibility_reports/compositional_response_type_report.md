# Compositional / Multinomial-Share Response Type Report

## Scope

This report evaluates how difficult it would be to add a new
`response_type = "compositional"` to the package, where "compositional" means:

- a vector of `K >= 3` non-negative shares per subject that sum to a fixed
  total (conventionally 1)
- examples: vote-share across `K` parties, market-share allocation across
  `K` products, budget/time-allocation across `K` categories, compositional
  chemistry assay outputs

This is different from the existing scalar `proportion` response type, which
covers exactly one bounded-`[0,1]` outcome per subject. Confirmed by reading
every concrete `proportion` inference class:
`EDI/R/inference_proportion_beta.R:33`,
`EDI/R/inference_proportion_fractional_logit.R:33`,
`EDI/R/inference_proportion_gcomp.R:84`, and
`EDI/R/inference_proportion_zero_one_inflated_beta.R:42` all call
`assertResponseType(des_obj$get_response_type(), "proportion")` and every one
fits a single scalar response — there is no existing vector-response
`proportion` path to extend.

This report is the compositional-specific companion to
`package_metadata/multivariate_response_type_report.md` (general
vector-valued outcomes, drafted separately). Compositional data **is** a
special case of multivariate data, but the simplex constraint
(`sum_k y_k = 1`, `y_k >= 0`) gives it extra structure that a well-known
transform can exploit — which is the central finding of this report and the
reason its difficulty verdict differs from the fully general multivariate
case.

The report covers:

- how common compositional outcomes are in experimental literatures
- the core package plumbing needed to admit a new response type that is a
  fixed-length vector rather than a scalar
- whether a log-ratio transform lets compositional outcomes reuse EDI's
  existing scalar-continuous machinery almost unchanged, or whether genuinely
  new multivariate/simplex-aware C++ kernels (Dirichlet regression) are
  required
- a pragmatic set of compositional-specific inference paths to implement first

## Short Answer

Adding the **label** `compositional` is easy, exactly as with `nominal`.

Adding a **coherent compositional response type** is **hard**, but
importantly *not uniformly hard* — the difficulty splits cleanly along one
axis:

- **The log-ratio-transform route (ILR/ALR/CLR) is moderate.** Once the
  simplex vector is transformed into `K-1` unconstrained real coordinates,
  each coordinate is an ordinary `continuous` response and can, in principle,
  reuse EDI's *existing* `InferenceContinOLS` / KK-continuous fitting
  machinery almost verbatim per-coordinate. The genuinely new work is mostly
  plumbing: vector-valued response storage, the transform/inverse-transform
  layer, and a joint (not per-coordinate-only) test/estimand.
- **The Dirichlet-regression route is hard.** It requires a new,
  simplex-native C++ MLE kernel (multivariate digamma/trigamma machinery,
  a `(K-1)`-dimensional score and Hessian, a genuinely new likelihood-test
  spec shape), which is comparable in scope to adding any other brand-new
  likelihood family to the package (see e.g.
  `package_metadata/quantile_regression_cpp_kernel_spec.md` for what a
  from-scratch kernel spec of similar size looks like).

Both routes require the same foundational architecture change: EDI's
response storage is a single scalar `numeric` vector on every design object
(`private$y`, confirmed at `EDI/R/design_abstract.R:124`,
`EDI/R/design_abstract.R:226`, `EDI/R/design_abstract.R:426`) — there is no
existing multi-column response storage anywhere in the package. That
foundational change is shared work, done once, that both routes build on.

**Recommended path:** implement the log-ratio-transform route first
(moderate difficulty, reuses existing scalar machinery, ships value quickly),
and treat true Dirichlet regression as a possible hard follow-on, not a
prerequisite.

## How Common Are Compositional Outcomes In Experimental Literatures?

Per `package_metadata/references/response_types_landscape_report.md`'s
"Additional Response Families Worth Noting" section (`## 5. Compositional /
multinomial-share outcomes`), the motivating examples are:

- market-share allocations
- vote-share vectors across parties (political-experiments field section,
  `response_types_landscape_report.md:361-392` — vote choice / party
  affiliation appear there as a `nominal` example, but the *share* version,
  i.e. the aggregate vote-share vector across a set of candidates rather than
  one voter's discrete choice, is the compositional analog)
- compositional chemistry outputs
- budget-share or time-allocation outcomes

### Political science and public opinion

Compositional analysis is a recognized subfield of political methodology —
aggregate vote-share data (by precinct, district, or survey cross-tab) is a
canonical compositional-data example in the broader statistics literature
(Aitchison's foundational compositional-data framework was popularized partly
through such applications). Field and lab experiments that report *aggregate*
treatment effects on vote or preference share across more than two options
(rather than one subject's single discrete choice) are the direct match for
this response type.

### Economics and marketing

Market-share response to advertising/pricing treatments, and budget or
time-allocation outcomes (e.g. how subjects reallocate a fixed budget across
`K` goods after a treatment), are natural compositional outcomes in applied
micro and marketing-science experiments.

### Chemistry / materials science

Compositional chemistry assay outputs (mixture/blend proportions) are a
classical compositional-data application outside the social sciences,
supporting the "present but specialized, not dominant" framing used
throughout the landscape report.

### Practical interpretation for this package

Compositional outcomes are a real but narrow slice of what EDI's target
audience would want — narrower than `nominal`, which arises any time a
per-subject discrete outcome is unordered. Compositional outcomes require the
outcome *itself* to be a vector of shares, which typically only arises when
the unit of observation is an aggregate (a district, a firm, a time period) or
when a single subject is asked to allocate a fixed resource across several
categories. This narrows the audience relative to `nominal` but does not make
it niche: political-science and marketing experiments both plausibly want it.

## What Already Exists

Grepping the entire proportion family and the generic dispatch layer
confirms there is currently:

- **No vector-valued response storage.** `Design`'s response field is a
  single `numeric` vector, one value per subject
  (`EDI/R/design_abstract.R:124,226,426`). There is no `y_matrix` or
  equivalent anywhere in `EDI/R/design_abstract.R`.
- **No simplex/compositional response_type value.** `response_type` dispatch
  (`EDI/R/globals.R`, and every `assertResponseType()` call site) recognizes
  exactly six values: `continuous`, `incidence`, `count`, `proportion`,
  `survival`, `ordinal` — confirmed identically in
  `package_metadata/references/response_types_landscape_report.md`'s own
  "Implementation TODOs" verification. `compositional`/`nominal` are absent
  from every dispatcher.
- **No log-ratio transform utilities.** There is no ILR/ALR/CLR helper
  anywhere in `EDI/R/*.R` (grepped for `ilr`, `alr`, `clr`, `compositional`,
  `simplex` — no hits outside this new report).
- **A close structural analog exists in `transform_cont_y_based_on_response_type()`**
  (`EDI/R/simulations_framework.R:126-150`), which already converts a raw
  continuous latent draw into a `response_type`-specific scale for
  simulation purposes. This is the right place to add a compositional branch
  later (generate `K` correlated continuous latents, then map through the
  inverse-ILR onto the simplex) — but it currently only ever emits a single
  scalar column per subject, confirmed by its `switch(response_type, ...)`
  structure over the six known types.
- **`InferenceMixinKKGeeShared`** (`EDI/R/inference_mixin_kk_gee_shared.R`)
  already threads a `gee_response_type` distinct from `des_obj$get_response_type()`
  in places (`inference_mixin_kk_gee_shared.R:78,143`), showing the codebase
  has *some* precedent for an inference class privately handling multiple
  outcome "shapes" under one banner — a useful but limited precedent, since
  it is still ultimately fitting one scalar GEE response per call, not a
  vector.

## Difficulty By Layer

### 1. Design base classes: moderate

`Design`'s response storage (`EDI/R/design_abstract.R:124,226,426`) needs a
second code path for vector responses: either (a) store `y` as an
`n x K` matrix instead of a length-`n` vector when `response_type ==
"compositional"`, or (b) store `K` separate scalar `y` columns and thread a
`K` count through the class. Option (a) is architecturally cleaner but
touches every method in `design_abstract.R` that currently assumes
`length(private$y) == n` (subsetting for jackknife/bootstrap replicate
construction, `private$y = private$y[i_b]` at
`EDI/R/design_abstract.R:490`, would need to become row-subsetting). This is
mechanical but not small — every `Design*` subclass that reads or writes `y`
needs an audit.

### 2. Most design classes: easy

The randomization/allocation logic in concrete `Design*` subclasses
(Bernoulli, blocking, IBCRD, KK14/KK21, etc.) does not depend on the response
at design time — response values are attached after allocation. So once the
base-class storage layer accepts a matrix response, individual design
subclasses need no changes, mirroring the same finding in
`package_metadata/nominal_response_type_report.md`'s "Most design classes:
easy" section.

### 3. KK21 / response-adaptive weighting designs: moderate to hard

KK21-style designs (`EDI/R/design_seq_one_by_one_KK21.R`,
`design_seq_one_by_one_KK21_stepwise.R`) use the *accruing response* to
influence future allocation. A compositional response has no single natural
scalar summary to feed a response-adaptive rule (unlike `incidence`/`count`,
which map directly onto a scalar imbalance or an outcome-adaptive
randomization statistic). Any KK21-with-compositional-outcome design would
need an explicit user-chosen scalarization (e.g. imbalance in one focal
component, or a distance-to-target-composition metric) — this is a genuine
open design decision, not just an implementation task, and should be
explicitly out-of-scope for a first pass (mirrors the same "moderate to hard"
verdict and reasoning as the nominal report's KK21 section).

### 4. `SimulationFramework`: hard

`SimulationFramework`'s core scalar assumptions run deep:
`transform_cont_y_based_on_response_type()` emits one column
(`EDI/R/simulations_framework.R:126-150`), and downstream truth/summary
metrics are built around a single `beta_T` scalar treatment effect (MSE and
coverage are both defined against a scalar `beta_T`,
`EDI/R/simulations_framework.R:204-206`). Compositional data generation needs
`K` correlated latents mapped through an inverse-ILR, and the ground-truth
"treatment effect" itself is ambiguous for a vector outcome (is it a
per-component effect vector? a single log-ratio contrast? a distance metric
between treatment and control mean compositions?) — this is the same open
estimand question flagged in "The Main Conceptual Problem" below, and it
must be resolved before `SimulationFramework` support is meaningful, not just
before it is implemented.

## The Main Conceptual Problem: What Is The Compositional Estimand?

This is the crux, exactly as "What Is The Nominal Estimand?" is the crux of
the nominal report. A compositional outcome does not have one natural
treatment-effect estimand; candidates include:

1. **Per-component contrasts**, e.g. `E[y_k | T=1] - E[y_k | T=0]` for a
   scientifically focal component `k` (e.g. the incumbent's vote share).
   Simplest to interpret, but ignores the sum-to-one constraint (a change in
   `y_k` mechanically implies changes elsewhere) and does not naturally
   support Type-I-error control across `K` components without an explicit
   multiplicity strategy.
2. **A log-ratio contrast**, e.g. the treatment effect on
   `ilr(y)` or a specific `alr(y)_k = log(y_k / y_ref)` coordinate. This is
   the natural compositional-statistics answer (Aitchison geometry treats the
   simplex as its own vector space where log-ratios are the "linear"
   coordinates) and is exactly what lets the log-ratio route reuse EDI's
   *existing* scalar-continuous machinery per transformed coordinate — but it
   is a coordinate-dependent, less intuitive quantity to report to an
   applied user than "the treatment moved vote share by X points."
3. **A global distributional/distance test**, e.g. whether the mean
   composition under treatment differs from control at all (an
   ILR-transformed Hotelling's `T^2` or permutation test on the full vector),
   with no single scalar effect size — structurally the compositional analog
   of `InferenceAllSimpleWilcox`'s or a global nominal test's "no natural
   single-number effect" problem, discussed in
   `package_metadata/nominal_response_type_report.md:476-501`
   ("Global two-sample nominal-distribution test").
4. **Dirichlet-regression coefficients**, i.e. `K-1` regression coefficient
   vectors (one per non-reference component) from a genuine simplex-native
   MLE, each interpretable as a log-odds-ratio-style effect relative to the
   reference component — the most statistically principled option and the
   one requiring a brand-new C++ kernel.

Unlike the nominal report, where the report ultimately recommends starting
with per-category contrasts plus a global test, compositional data has a
cleaner default answer from the statistics literature: **the log-ratio
route (option 2) is the standard approach** for exactly this reason — it
converts an otherwise ambiguous vector estimand into `K-1` ordinary scalar
estimands with a well-defined inverse transform back to interpretable
compositional space. This report recommends leading with option 2, offering
option 1 (per-component contrasts, with an explicit multiplicity-adjustment
flag) as a simpler complementary output, and treating options 3 and 4 as
later, harder additions.

## How The Existing `inference_all_*` Paths Would Respond

- **`InferenceAllSimpleMeanDiff`** (mean-difference base class): would need
  to reject a raw `compositional` response outright, exactly as it rejects
  `nominal` today (`package_metadata/nominal_response_type_report.md:297-317`)
  — a vector has no scalar mean difference. It becomes directly usable,
  unmodified, on each ILR-transformed coordinate once the transform layer
  exists, since each coordinate is an ordinary continuous scalar.
- **`InferenceAllSimpleWilcox`**: same story — rejects the raw vector,
  becomes usable per-transformed-coordinate.
- **`InferenceMLEorKMSummaryTable`**: this abstract layer is response-shape
  agnostic per the nominal report's own finding
  (`nominal_response_type_report.md:333-350`) and should require no changes
  to remain neutral for `compositional` either.
- **`InferenceAsympLikStdModCache`**: directly reusable, unmodified, for
  fitting each ILR-transformed coordinate as an ordinary continuous
  likelihood model (this is the single biggest reuse win of the log-ratio
  route: `K-1` calls into machinery that already exists and is already
  optimized). It is *not* reusable as-is for true Dirichlet regression, which
  needs a joint `(K-1)`-parameter-block likelihood, score, and Hessian that
  this cache layer does not currently know how to shape.
- **Randomization / bootstrap base classes**: mostly neutral, mirroring the
  nominal report's finding — resampling operates on subject indices, not on
  the response's shape, so a vector response subsets the same way a scalar
  one does once the storage layer supports it (`design_abstract.R:490`'s
  `private$y = private$y[i_b]` becomes a row-subset instead of an
  element-subset).
- **`InferenceSuite`**: will need to explicitly branch on `compositional`
  wherever it currently assumes a single scalar summary is being reported
  per fitted path (table columns, plotting) — an audit-and-guard pass, not a
  redesign, once the response-shape and estimand decisions above are locked.

## Package Areas That Need Explicit `compositional` Fences

Mirroring the nominal report's approach, every one of these needs either an
explicit `compositional` branch or an explicit rejection with a clear error,
so a vector response never silently gets treated as `NA`/coerced/truncated
by code that assumes a scalar:

- `Design` response storage and subsetting (`design_abstract.R`)
- `assertResponseType()` and the six-way `response_type` enum
  (`EDI/R/globals.R`)
- `transform_cont_y_based_on_response_type()`
  (`EDI/R/simulations_framework.R:126-150`)
- `SimulationFramework`'s truth/summary metric machinery (scalar `beta_T`
  assumption, `simulations_framework.R:204-206`)
- Every generic "all-subject" `inference_all_*` base class listed above
- `InferenceSuite`'s reporting/plotting layer

## Common Compositional-Specific Inference Paths To Implement

### 1. ILR-transformed per-coordinate continuous regression (the core deliverable)

Apply an isometric log-ratio transform to map the `K`-simplex response onto
`K-1` unconstrained real coordinates, fit each coordinate with EDI's existing
`InferenceContinOLS` (or another existing continuous path) machinery, and
report both the per-coordinate estimates/CIs and their back-transform into
an interpretable compositional effect (e.g. the implied shift in each
component's expected share holding the others' ratios fixed). This is the
single highest-leverage path: nearly all its numerical core is *already
implemented and already fast* (`EDI/R/inference_continuous_ols.R`), so the
new work is almost entirely the transform layer, vector storage, and a
combining/reporting wrapper — not new numerical kernels.

### 2. Global joint test on the transformed vector

An ILR-transformed Hotelling's `T^2` (or a permutation-based multivariate
generalization, reusing the package's existing randomization-inference
scaffolding) testing whether the mean composition differs between arms at
all, without requiring a scalar effect size — the compositional analog of
the nominal report's "Global two-sample nominal-distribution test".

### 3. Per-component contrast with multiplicity adjustment

Report simple `E[y_k|T=1] - E[y_k|T=0]` contrasts for each of the `K`
components with a Holm/Bonferroni-style adjustment flag, as a
simpler-to-interpret complement to the ILR-based estimand for applied users
who want "did component `k`'s share change" rather than a log-ratio
coordinate. Check whether any existing multiplicity-adjustment utility
already exists in the package before writing a new one (grepped
`EDI/R/*.R` for `holm`/`bonferroni`/`p.adjust` — none found, so this would be
new, small utility code, not a numerical kernel).

### 4. Dirichlet regression (later, hard follow-on)

A genuine `K`-component-simplex MLE via Dirichlet regression, requiring a
new C++ kernel (multivariate digamma/trigamma-based score and Hessian,
comparable in scope to `package_metadata/quantile_regression_cpp_kernel_spec.md`).
Recommended as an explicit later phase, not part of the first wave.

## Recommended First-Wave Compositional Inference Set

1. ILR-transformed per-coordinate continuous regression (path 1) — ships the
   most value for the least new numerical work.
2. Global ILR-transformed joint test (path 2) — cheap to add once the
   transform layer and per-coordinate fits exist, and answers the "is there
   any effect at all" question applied users will ask first.
3. Per-component contrasts with multiplicity adjustment (path 3) — small,
   complements 1 and 2 with a more intuitive report format.

Defer Dirichlet regression (path 4) to a later phase explicitly informed by
user demand, given its cost is comparable to a full new likelihood family.

## `SimulationFramework` Implications

### Data generation

Add a `compositional` branch to `transform_cont_y_based_on_response_type()`
(`simulations_framework.R:126-150`) that generates `K` correlated Gaussian
latents (reusing whatever correlated-draw machinery the framework already
has for other multi-parameter response types, if any exists — otherwise a
new small helper) and maps them through the inverse-ILR onto the simplex,
parameterized by a baseline composition and a treatment-effect shift applied
in log-ratio space.

### Treatment effect semantics

The framework's `beta_T` must become either (a) a `K-1`-length vector of
per-log-ratio-coordinate effects (natural fit for the recommended estimand in
"The Main Conceptual Problem" above, but breaks every scalar-`beta_T`
assumption in the MSE/coverage machinery at
`simulations_framework.R:204-206`), or (b) restricted, for a first pass, to a
single scalar shift applied uniformly to one designated log-ratio coordinate
(much smaller change, ships faster, sacrifices generality). Recommend (b)
for the first wave, matching the "ship the log-ratio route first" theme of
this whole report.

### Truth and summary metrics

MSE and coverage need a `compositional`-aware definition once `beta_T` is
vector-valued (b) above defers this problem by keeping `beta_T` scalar for
the first pass, so no change is needed there initially; a later expansion to
the full `K-1`-vector case (a) would need genuinely new per-coordinate or
Euclidean/Aitchison-distance-based summary metrics.

## Recommended Implementation Plan

### Stage 1: Admit the type and vector storage safely

Add `compositional` to the `response_type` enum and `assertResponseType()`
dispatch (`EDI/R/globals.R`), extend `Design`'s response storage to accept an
`n x K` matrix when `response_type == "compositional"`
(`EDI/R/design_abstract.R`), and audit every subsetting/replicate-generation
site that currently assumes `length(private$y) == n`. Every other
`inference_all_*` base class should explicitly reject `compositional` at
this stage (fence, don't silently mishandle).

### Stage 2: Build the ILR transform layer

Add ILR/inverse-ILR helper functions (new, small, no external heavyweight
dependency needed — the transform is closed-form linear algebra). Wire them
into a new `InferenceCompositionalILR`-style wrapper class that: transforms
the stored `n x K` response into `n x (K-1)` continuous coordinates,
delegates each coordinate's fit to the existing `InferenceContinOLS` (or
another chosen continuous path), and inverse-transforms results back to
compositional space for reporting.

### Stage 3: Add the global joint test

Layer a Hotelling's-`T^2`-style (or permutation-based) global test on top of
the Stage 2 per-coordinate fits, reusing existing randomization-inference
scaffolding where possible.

### Stage 4: Add per-component contrasts with multiplicity adjustment

Small, additive; a simpler complementary report format for applied users.

### Stage 5: Extend `SimulationFramework`

Add the `compositional` data-generation branch (scalar-`beta_T`-restricted
first pass per the "Treatment effect semantics" discussion above), then test
the Stage 2-4 inference paths against it end-to-end.

### Stage 6 (later, separately scoped): Dirichlet regression

Only after real user demand is established; treat as its own project on the
scale of `package_metadata/quantile_regression_cpp_kernel_spec.md`, not a
continuation of Stages 1-5.

## Bottom Line

Adding `compositional` support to EDI is **hard as a fully general feature**
(true Dirichlet regression, vector `beta_T`, all four inference paths) but
**moderate for the recommended first wave** (ILR transform + reuse of
existing scalar-continuous machinery + a global test + per-component
contrasts). The foundational architecture change — vector-valued response
storage on `Design` — is required either way and is itself a moderate,
mechanical (not conceptual) piece of work. Unlike `nominal`, where the
estimand question has several genuinely competing candidates with no clean
default, compositional data has a well-established statistical answer (the
log-ratio transform) that converts most of the problem back into work EDI's
existing continuous-response machinery already does well. That is the
central reason this report's verdict is more optimistic than the raw `hard`
label in the landscape survey would suggest on its own — the `hard` label is
accurate for the fully general feature, but the recommended first wave is
substantially cheaper than that label implies.
