# Rank / Preference / Choice-Set Response Type Report

## Scope

This report evaluates how difficult it would be to add first-class support for
**rank / preference / choice-set outcomes** to the package. This covers two
related but statistically distinct outcome shapes:

- **(a) single discrete choice**: a subject selects one option from a
  candidate/option set (e.g. a discrete-choice experiment, or "which
  candidate would you vote for")
- **(b) full or partial rankings**: a subject orders some or all of the
  options from most to least preferred (e.g. ranked-choice ballots,
  best-worst scaling, product ranking tasks)

This report fulfills TODO-6 in
[response_types_landscape_report.md](references/response_types_landscape_report.md),
which flagged this family as `hard`, motivated by the "Political experiments"
field section's ranked-candidate and discrete-choice examples (see
[response_types_landscape_report.md:641-662](references/response_types_landscape_report.md:641)).
It does not cover ordinal preference *ratings* (e.g. a 1-5 Likert scale) —
that is squarely the existing `ordinal` response type — only outcomes that are
a *choice among*, or an *ordering of*, a set of discrete alternatives.

The central finding of this report is that **(a) and (b) are not the same
project**. Single discrete choice is architecturally almost identical to the
`nominal` response type already assessed in
[nominal_response_type_report.md](nominal_response_type_report.md); full/partial
rankings are a genuinely different statistical object that the package's
current likelihood-test architecture cannot host without new machinery. The
two sub-cases get separate verdicts below.

## Short Answer

- **(a) Single discrete choice**: **moderate**, effectively riding on the
  `nominal` response type's plan. If `nominal` is added first (per
  [nominal_response_type_report.md](nominal_response_type_report.md)), discrete
  choice needs almost no extra architecture — it *is* a nominal outcome. If
  `nominal` is not added first, this report's Stage 1/2/3 below duplicate that
  work.
- **(b) Full/partial rankings**: **hard**. The outcome is not a single
  category draw but an ordered permutation (or partial permutation) of the
  option set. This breaks the package's `Design`-level assumption of one
  scalar-or-factor response per subject, and breaks the likelihood-test
  architecture's assumption that a fitted model produces a coefficient vector
  `b` with a designated scalar treatment entry `b[2]` (see "The Main
  Conceptual Problem" below). A Plackett-Luce or Mallows model is a new
  numerical core, not a new `response_type` branch.

## How Common Are Rank/Choice Outcomes In Experimental Literatures?

[response_types_landscape_report.md:641-662](references/response_types_landscape_report.md:641)
already summarizes this: "common enough in political science, economics, and
marketing-style experiments to matter, but too structurally different to be
folded cleanly into the current scalar response menu." This report does not
re-derive that literature survey; it inherits the verdict and focuses on the
architectural question of *how* structurally different the two sub-cases are
from each other and from the existing package.

Concretely, per the landscape report's "Political experiments" field section
([response_types_landscape_report.md](references/response_types_landscape_report.md)),
vote choice and party affiliation are already-cited motivating examples for
`nominal` — those are **single discrete choice** cases. Ranked-candidate
preference (e.g. ranked-choice voting, conjoint/discrete-choice experiments
with rank-ordered responses) is the **(b)** case.

## What Already Exists

The package has no response type for either sub-case today. The closest
existing analogs are:

- `ordinal`: an *ordered* categorical outcome with a single natural scale
  direction — see [design_abstract.R:94](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_abstract.R:94),
  which lists the full allowed `response_type` set as `"continuous",
  "incidence", "proportion", "count", "survival", "ordinal"` — note neither
  `nominal` nor any rank/choice type is present, confirming both are entirely
  new additions today.
- `InferenceOrdinalAdjCatLogit` / `InferenceOrdinalStereotypeLogit`
  ([inference_ordinal_adj_cat_logit.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_ordinal_adj_cat_logit.R),
  [inference_ordinal_stereotype_logit.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_ordinal_stereotype_logit.R)):
  these already fit multi-parameter, per-category logit models and are the
  closest existing template for *any* new choice-among-categories model,
  because they demonstrate the package already knows how to fit models with
  more free parameters than one scalar treatment effect and reduce the output
  to a single `beta_hat_T` for the package's asymptotic CI/p-value API (see
  `get_likelihood_test_spec()` at
  [inference_ordinal_stereotype_logit.R:170](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_ordinal_stereotype_logit.R:170)).
- Nothing in `EDI/src/*.cpp` implements a multinomial-logit, conditional-logit,
  or rank-ordered-logit likelihood today (confirmed via
  `grep -rli "multinom\|plackett\|rank.*logit\|conditional.*logit" EDI/src/*.cpp`
  returning no matches other than the KK21 weighting kernels' unrelated
  "logistic" naming).

## Sub-Case (a): Single Discrete Choice

### Relationship to `nominal`

A single discrete choice among `K` unordered alternatives, observed once per
subject, is definitionally a nominal outcome. There is no additional
structure beyond what
[nominal_response_type_report.md](nominal_response_type_report.md) already
analyzes in full. Specifically, from that report:

- **Difficulty by layer** (design base classes: moderate; most design classes:
  easy; KK21 response-adaptive weighting: moderate to hard; `SimulationFramework`:
  hard) — see
  [nominal_response_type_report.md:176-261](nominal_response_type_report.md:176)
  — carries over **verbatim**. Nothing about "this nominal outcome happens to
  represent a choice among options" changes any of those four layers.
- **The main conceptual problem** — "what is the nominal estimand?" — carries
  over verbatim (
  [nominal_response_type_report.md:263-293](nominal_response_type_report.md:263)).
  A discrete-choice experiment adds one simplifying detail: the choice set is
  often the *design's own treatment arms or a small enumerated option list*,
  which makes "focal-category vs. rest" or "multinomial logit with a natural
  baseline option" an even more natural first estimand than in the generic
  nominal case.
- **Recommended first-wave set**
  (`InferenceNominalPearsonChisq`, `InferenceNominalCategoryContrast`,
  `InferenceNominalMultinomLogit`,
  [nominal_response_type_report.md:520-544](nominal_response_type_report.md:520))
  directly answers discrete-choice questions: "did treatment change the
  probability of choosing option k" is exactly
  `InferenceNominalCategoryContrast`, and "did treatment shift the full choice
  distribution" is exactly `InferenceNominalMultinomLogit`.

### What discrete choice adds beyond generic nominal

Two things are genuinely choice-specific and not covered by the nominal
report, both minor:

1. **Choice-set size/composition varies across subjects (menu variation)** —
   if different subjects face different option menus (a standard conjoint /
   discrete-choice-experiment design), the per-subject choice probability
   model becomes a *conditional* logit (probabilities normalized over the
   subject's own offered menu) rather than a *multinomial* logit
   (probabilities normalized over a fixed global category set). This is a
   real difference in the C++ likelihood kernel (conditional logit needs a
   per-subject menu/offer-set input, not just a global category label) but is
   still a scalar-treatment-effect, single-choice-per-subject likelihood —
   it does not disturb `Design`, `SimulationFramework`, or the asymptotic-CI
   architecture any differently than plain multinomial logit does.
2. **No new package layer is needed if choice sets are fixed** (every subject
   sees the same K options) — this degenerates exactly to nominal
   `InferenceNominalMultinomLogit`.

### Verdict: (a) Single discrete choice — moderate

Same verdict as generic `nominal`. If `nominal` is implemented first, discrete
choice is a documentation/labeling exercise plus (optionally) a conditional-
vs-multinomial kernel variant for menu-varying designs. If implemented
standalone, it requires redoing Stages 1-4 of the nominal report's plan.

## Sub-Case (b): Full / Partial Rankings

### The Main Conceptual Problem: What Is A Ranking, Architecturally?

Every response type currently in the package — including the newly-analyzed
`nominal` case — is stored as **one scalar or one factor level per subject**
in the `Design` base class's `y` vector
(`private$y = rep(NA_real_, n)`,
[design_abstract.R:123](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_abstract.R:123)),
and every value-add/assert method (`add_one_subject_response()`,
`assert_y()`, etc.) is written assuming a subject contributes exactly one
number or one factor level.

A ranking is not a scalar or a factor level — it is a **permutation (or
partial permutation) over a fixed option set**, i.e. one integer vector of
length up to `K` per subject (or, for partial rankings, a shorter ordered
subset). This means:

- the `Design$y` storage contract itself (one `NA_real_`/factor cell per
  subject) cannot represent a ranking without a redesign of the response
  container — this is a strictly larger change than anything the nominal
  report proposed, because nominal only needed a new *factor* semantics, not
  a new *response shape*.
- `assert_y()`, `add_one_subject_response()`, and
  `add_all_subject_responses()` (all response-shape-aware entry points
  enumerated for nominal in
  [nominal_response_type_report.md:186-193](nominal_response_type_report.md:186))
  would all need a genuinely new code path, not an extra `switch` arm, because
  their current signatures accept one value/row, not a variable-length ordered
  set.

### Likelihood-Test Architecture: Does It Fit?

The package's shared likelihood-backed inference layer,
`inference_asymp_lik_std_mod_cache` in
[inference_all_abstract_asymp_lik_std_mod_cache.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_asymp_lik_std_mod_cache.R),
assumes:

- `generate_mod()` returns a fitted-model object from which a scalar treatment
  estimate `beta_hat_T` and its standard error are extracted (see
  `compute_estimate()` at
  [inference_all_abstract_asymp_lik_std_mod_cache.R:9](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_asymp_lik_std_mod_cache.R:9)
  and `get_standard_error()` at
  [inference_all_abstract_asymp_lik_std_mod_cache.R:52](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_asymp_lik_std_mod_cache.R:52)),
  used by 50 concrete `get_likelihood_test_spec()` implementations across the
  package (`grep -rln "get_likelihood_test_spec" EDI/R/*.R | wc -l` → 50).

A Plackett-Luce or Mallows-type rank model **can** in principle be fit by
maximum likelihood and **can** in principle expose one scalar treatment
coefficient (e.g. "treatment shifts the log-worth of the treated arm's
option by `beta_T`"), so this layer is not structurally incompatible the way
`Design$y` storage is. The blocker is not `get_likelihood_test_spec()` itself
— it is that no C++ numerical core exists to compute a Plackett-Luce (or
Mallows) log-likelihood, score, and Hessian, and building one is a new
optimization project of comparable scope to (say) the existing ordinal
adjacent-category/stereotype logit kernels, not a thin wrapper. Concretely,
the Plackett-Luce log-likelihood for a full ranking of `K` items is a sum of
`K-1` sequential multinomial-choice terms conditioned on the remaining
unranked items — this needs the same conditional-logit kernel called out in
sub-case (a)'s "menu variation" point, applied recursively per rank position,
which compounds the computational cost per subject roughly `K`-fold relative
to a single discrete choice.

### `SimulationFramework`: harder than nominal

`SimulationFramework`'s `transform_cont_y_based_on_response_type()`
([simulations_framework.R:126-152](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R:126))
transforms one latent continuous draw `y_cont` per subject into one observed
response value via a `switch(response_type, ...)`. Even the hardest existing
case (`ordinal`, via quantile-cut) still produces one scalar per subject. A
ranking generator needs `K` latent utilities per subject (one per option,
e.g. via a Plackett-Luce generative model: draw independent Gumbel/exponential
noise per option, rank by realized utility) and returns a length-`K`
permutation, not a scalar — this is a strictly different generator shape, not
a new `switch` arm, and it also breaks the framework's scalar `betaT`
semantics and scalar truth/MSE/coverage/power summary machinery
(same problem the nominal report identifies at
[nominal_response_type_report.md:578-593](nominal_response_type_report.md:578),
but without nominal's fallback of "target one focal category's log-odds" —
for a full ranking there is no obviously-scalar analog short of picking a
single focal-item's worth parameter, which is a reasonable simplification but
throws away most of the ranking's information).

### Verdict: (b) Full/partial rankings — hard

This is genuinely new statistical territory for the package: a new response
container shape at the `Design` layer, a new C++ likelihood kernel family
(Plackett-Luce / conditional-logit-recursive), and a new `SimulationFramework`
generator. It is comparable in scope to adding an entirely new likelihood
family (e.g. on the order of the existing ordinal-regression kernel suite),
not a `response_type` label addition.

## Package Areas That Need Explicit Fences

If discrete choice (a) is implemented via `nominal`, the same fencing
requirements apply as in
[nominal_response_type_report.md:395-414](nominal_response_type_report.md:395)
(`InferenceAllSimpleMeanDiff`, `InferenceAllSimpleMeanDiffPooledVar`,
`InferenceAllSimpleWilcox`, quantile-regression abstractions, and any
response-type-agnostic path lacking `assertResponseType(...)`).

If full/partial rankings (b) are ever implemented, they would need a **new**
`response_type` value (e.g. `"ranking"`) distinct from `nominal`, because a
ranking is not representable in the same `y` storage — this means every
currently-nominal-tolerant guard would need to separately consider and likely
reject `"ranking"` as well, doubling the fencing surface relative to nominal
alone.

## Recommended Implementation Plan

### Stage 1 (shared with `nominal`): Admit discrete-choice as `nominal`

Implement `nominal` per
[nominal_response_type_report.md](nominal_response_type_report.md)'s Stages
1-3. This alone delivers useful single-discrete-choice inference
(`InferenceNominalCategoryContrast` answers "did treatment change probability
of choosing option k", `InferenceNominalMultinomLogit` answers the full
choice-distribution question).

### Stage 2: Conditional-logit variant for menu-varying discrete choice

Only needed if the package wants to support designs where subjects see
different option menus (standard in conjoint/discrete-choice experiments).
Add a conditional-logit C++ kernel alongside the multinomial-logit kernel from
Stage 1, taking a per-subject offer-set indicator as an additional input.

### Stage 3: Ranking response type (only if full/partial rankings are a real
package goal)

This is a standalone project, not a Stage 4 of the nominal or discrete-choice
work:

1. Design a new `y`-storage shape for permutation-valued responses at the
   `Design` layer (this is the largest single piece of new plumbing).
2. Implement a Plackett-Luce (or Mallows) log-likelihood, score, and
   information/Hessian C++ kernel, exposing one scalar focal-item treatment
   coefficient for `get_likelihood_test_spec()` compatibility.
3. Add a `SimulationFramework` ranking generator (per-option latent utilities
   → realized ranking) and a scalar truth/`betaT` definition (e.g. treatment
   effect on the focal item's Plackett-Luce worth parameter).
4. Add explicit rejection fences to every generic `inference_all_*` path that
   would otherwise silently misinterpret a ranking's storage as ordinary
   numeric/factor data.

## Bottom Line

- **Single discrete choice**: **moderate** — this is the `nominal` response
  type wearing a different hat. Implement `nominal` first
  ([nominal_response_type_report.md](nominal_response_type_report.md)) and
  discrete choice comes almost for free; only menu-varying designs need an
  extra conditional-logit kernel (Stage 2).
- **Full/partial rankings**: **hard** — a new response-storage shape at the
  `Design` layer, a new Plackett-Luce/Mallows C++ likelihood core (recursive
  conditional-logit, no existing analog in `EDI/src`), and a new
  `SimulationFramework` generator with no scalar-truth fallback as clean as
  nominal's focal-category trick. Treat this as a standalone project on the
  scale of adding a new likelihood family, to be scoped only after Stage 1/2
  discrete-choice support (which reuses `nominal`) is in place and there is
  concrete demand for full rankings specifically.

The pragmatic recommendation is to **not** conflate the two sub-cases in
planning or in package API naming: ship discrete choice as part of `nominal`
work, and treat rank/partial-ranking support as a separate, later, harder
feasibility question — revisit only if applied demand for full rankings
(as opposed to single choice) actually materializes.
