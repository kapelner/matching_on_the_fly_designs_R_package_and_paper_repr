# IVWC Compound Class Consolidation — Design

## Background

The `Inference*` R6 hierarchy in `EDI/R` currently looks like (relevant subset):

```
InferenceAsymp
 -> InferenceMLEorKMSummaryTable   (Wald z/t inference off a coef+SE+df summary; no likelihood required at all)
     -> InferenceAsympLik          (adds wald/score/gradient/lik_ratio/bartlett dispatch via capability flags;
                                     spans full, partial, quasi, and composite likelihoods generically)
         -> InferenceParamBootstrap                    (coherent single likelihood, supports parametric-bootstrap LR)
             -> InferenceKKPassThroughCompound          (KK matched-pairs + reservoir compound, WITH param bootstrap)
         -> InferenceKKPassThroughCompoundNoParamBootstrap  (KK matched-pairs + reservoir compound, WITHOUT param bootstrap)
```

`InferenceKKPassThroughCompound` / `InferenceKKPassThroughCompoundNoParamBootstrap` (in
`inference_all_abstract_KK_passthrough_compound.R`) already centralize the logic needed
by "IVWC" (inverse-variance-weighted-combination) families: a matched-pairs estimate and
a reservoir estimate combined via a fixed weight `w_star`. The `NoParamBootstrap` variant
splices in `InferenceMixinKKPassThroughCompound`, whose private slot already sets
`supports_likelihood_tests = function() FALSE` and provides the shared
`compute_estimate_with_bootstrap_weights` fallback.

Eight files already build their IVWC classes on top of `InferenceKKPassThroughCompoundNoParamBootstrap`
correctly: `ols_ivwc`, `mean_diff_IVWC`, `quantile_rand_ci`, `bai_abstract`, `robust_regr_ivwc`,
`robust_regr_one_lik`, `newcombe_ivwc_univ`, `rank_regr_ivwc_abstract`.

## Problem

Thirteen other files inherit `InferenceAsympLik` directly instead, and hand-duplicate what the
shared compound base already provides:

- `private$supports_likelihood_tests = function() FALSE`
- `public$compute_asymp_confidence_interval` — `private$shared(); private$compute_z_or_t_ci_from_s_and_df(alpha)`
- `public$compute_asymp_two_sided_pval` — `private$shared(); private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)`
- `public$approximate_bootstrap_distribution_beta_hat_T = function(...) eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))`

Of those 13, six are confirmed genuine matched-pairs + reservoir compound estimators (they compute
`beta_m`/`beta_r` and combine via `w_star`) — i.e. they are exactly what
`InferenceKKPassThroughCompoundNoParamBootstrap` was built for, just not wired up to it. The other
three do not fit that shape at all (see Out of Scope below) and should not be forced onto it.

## Decision

No new class, and no rename of `InferenceMLEorKMSummaryTable`. (Both were considered and rejected —
see prior discussion: `InferenceAsympLik` already documents the full/partial/quasi/composite-likelihood
taxonomy, and `InferenceMLEorKMSummaryTable` correctly covers likelihood-free KM inference too, so
neither change would add a distinction that doesn't already exist.)

Instead: migrate the six confirmed-compound straggler classes to inherit from
`InferenceKKPassThroughCompoundNoParamBootstrap`, deleting the boilerplate they were duplicating.
This is a pure consolidation onto an existing, already-proven abstraction — not a new design.

The user has confirmed these six classes are legacy (kept only for theoretical comparison against
their superior `*OneLik` siblings), so exact-behavior preservation is not a hard constraint — but the
migration below is designed to be numerically inert anyway (same code, different inheritance path).

## In scope — migrate these six

| Class | File |
|---|---|
| `InferenceSurvivalKKStratCoxPHIVWC` | `inference_survival_KK_strat_cox.R` |
| `InferenceAbstractKKWeibullFrailtyIVWC` | `inference_survival_KK_weibull_frailty.R` |
| `InferenceSurvivalKKClaytonCopulaIVWC` | `inference_survival_KK_clayton_copula.R` |
| `InferenceIncidKKCondLogitIVWC` | `inference_incidence_KK_cond_logit.R` |
| `InferenceAbstractKKWilcoxBaseIVWC` | `inference_all_KK_wilcox_ivwc.R` |
| `InferenceAbstractKKLWACoxIVWC` | `inference_survival_KK_lwa_cox_ivwc_abstract.R` |

All six currently share the identical shape:

```r
InferenceXyzIVWC = R6::R6Class("InferenceXyzIVWC",
    lock_objects = FALSE,
    inherit = InferenceAsympLik,
    public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
        ...,
        compute_asymp_confidence_interval = function(alpha = 0.05){ ... },
        compute_asymp_two_sided_pval = function(delta = 0){ ... },
        approximate_bootstrap_distribution_beta_hat_T = function(...){ eval(body(...)) }
    ))),
    private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
        ...,
        supports_likelihood_tests = function() FALSE,
        ...
    )))
)
```

### Mechanical change, per file

1. `inherit = InferenceAsympLik` -> `inherit = InferenceKKPassThroughCompoundNoParamBootstrap`.
2. `public`: drop the `as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(...)))`
   wrapper entirely; keep a plain `list(...)` containing only the class-specific methods
   (`initialize`, `compute_estimate`, the `shared()`-driving logic, and any custom
   `compute_estimate_with_bootstrap_weights` override). Remove `compute_asymp_confidence_interval`,
   `compute_asymp_two_sided_pval`, and `approximate_bootstrap_distribution_beta_hat_T` — all now
   inherited from `InferenceKKPassThroughCompoundNoParamBootstrap`.
3. `private`: same de-wrapping; drop `supports_likelihood_tests = function() FALSE` (now inherited).
   Keep everything else (`compute_basic_match_data`, `shared()`, the `w_star` combination code, etc.).
4. If a class overrides `supports_lik_ratio_param_bootstrap()` (only `WeibullFrailtyIVWC` does,
   conditionally on `use_rcpp`) — leave it in place during migration but flag it for the follow-up
   cleanup below; it's provably inert on the `NoParamBootstrap` branch.

## Out of scope — leave inheriting `InferenceAsympLik` directly, unchanged

These do not have a matched-pairs + reservoir split and do not fit
`InferenceKKPassThroughCompoundNoParamBootstrap`'s shape:

- `InferenceOrdinalKKCondAdjCatLogitRegr` (`inference_ordinal_KK_cond_adj_cat_logit.R`) — single
  data-expansion conditional-logit fit over the whole sample, no pair/reservoir combination.
- `InferenceOrdinalPairedSignTest` (`inference_ordinal_paired_sign_test.R`) — matched-pairs-only
  design; bootstrap and jackknife are explicitly disallowed (design-constraint violation), so there
  is no reservoir component to combine.
- `InferenceAbstractKKOrdinalCLMM` (`inference_ordinal_KK_clmm_abstract.R`) — single joint CLMM fit
  over matched pairs and reservoir together (a mixed model, not a two-source weighted combination).

Also unchanged: `InferenceAsympLik`, `InferenceParamBootstrap`, `InferenceMLEorKMSummaryTable`,
`InferenceKKPassThroughCompound`, `InferenceKKPassThroughCompoundNoParamBootstrap` themselves — no
new methods needed on any of them.

## Expected behavior impact

None. This is a structural-only refactor: the same `shared()` implementations, the same
`w_star`-based combination arithmetic, and the same bootstrap-distribution body move from
"copy-pasted per file" to "inherited from one place." Numeric output (estimate, CI, p-value,
bootstrap distribution) for the six migrated classes should be bit-identical before and after.

## Testing plan

- Run the existing test suite (`comprehensive_tests.R` / targeted tests) for the six affected
  concrete classes before and after migration; diff estimate/CI/p-value output — expect no change.
- Re-run `Rscript fast_roxygenize.R` from the repo root afterward to regenerate NAMESPACE/docs,
  since the R6 inheritance chain is reflected in generated `.Rd` files (see e.g.
  `InferenceAbstractKKSurvivalRankRegrIVWC.Rd`'s inheritance diagram).
- `devtools::check()` / R CMD check to catch any method resolution issue from the flattened
  `inherit` chain.

## Considered and rejected

Three alternative restructurings came up while scoping this and were rejected:

1. **New `InferenceAsympPartialLik` class**, sitting above or below `InferenceAsympLik` to hold
   "partial-likelihood-only" methods. Rejected: `InferenceAsympLik` already documents itself as
   spanning full/partial/quasi/composite likelihoods generically via capability flags, and the
   actual duplication turned out to be the IVWC compound-estimator pattern addressed above, not a
   missing likelihood tier.

2. **Rename `InferenceMLEorKMSummaryTable` to something like `InferenceAsympFullPartialQuasiCompLik`**
   to describe the likelihood taxonomy it spans. Rejected: this class requires no likelihood at all
   (it's inherited by Kaplan-Meier-based inference, which is nonparametric) — its real contract is
   "produces a coef+SE+df summary supporting Wald z/t inference," which is broader and more accurate
   than any likelihood descriptor. The likelihood taxonomy is already correctly documented one level
   down, on `InferenceAsympLik`.

3. **Rename `InferenceParamBootstrap` to `InferenceFullLikParamBootstrap`**, to assert that every
   daughter implements a full likelihood. Rejected on a direct counterexample:
   `InferenceSurvivalKKStratCoxPHOneLik` (header: "Stratified Cox **Combined**-Likelihood Compound
   Inference," `shared_combined_likelihood()`) inherits `InferenceParamBootstrap` and sets
   `supports_lik_ratio_param_bootstrap = TRUE`, despite Cox's partial likelihood being the textbook
   example of a likelihood that is explicitly *not* full. "Full" is the wrong invariant; "a single
   likelihood coherent enough to simulate from" (i.e. `simulate_under_lik_null()` is implementable)
   is the real one, but that's cross-cutting rather than something a class name should assert.

4. **Restructure the whole hierarchy as a strict linear chain** `NoLikelihood < PartialLikelihood <
   FullLikelihood < ...`. Originally rejected here on the grounds that Cox PH (partial likelihood)
   reaches top capability, quasi-Poisson (quasi-likelihood) reaches zero capability despite being
   chained under `InferenceParamBootstrap`, GEE sits outside the whole subtree while taxonomically
   identical quasi-Poisson sits inside it, and `InferenceAbstractKKOrdinalCLMM` is a full joint
   likelihood flagged to zero capability — all read at the time as proof that taxonomy and capability
   are orthogonal, not nested.

   **⚠ Superseded.** On further pressure this conclusion doesn't hold: every one of those
   counterexamples conflates likelihood strength with two *different* things — whether an estimator
   is a compound of two separately-fit models (an orthogonal structural axis), and engineering gaps
   (an inconsistently-wired GEE, an unfinished R fallback path) that look like theory violations but
   aren't. Once those are separated out, a natural, textbook-grounded superset/subset chain does
   exist. See the new section below, which replaces this entry as the actual design going forward.

## Revised: a theory-grounded 4-tier likelihood taxonomy, with everything else as mixins

**Scope note:** everything below is a second, independent design, much larger in scope than the
six-file migration specified in "In scope" above (which remains valid, fully specified, and safe to
implement on its own). This section is the conceptual architecture only — it is not yet a
file-by-file implementation plan. It touches the abstract backbone and, transitively, every one of
the roughly 130 concrete inference kernels' `inherit =` line. See "Scope and sequencing" at the end
of this section for how the two relate.

### The four tiers

Grounded in the standard likelihood-methods spectrum (Cox 1975 partial likelihood; Wedderburn/
McCullagh quasi-likelihood; ordinary full-MLE theory — the same spectrum covered in e.g. Severini's
*Likelihood Methods in Statistics*). Each tier is a strict statistical and syntactic superset of the
one below it:

0. **`InferenceNoLik`** — no likelihood object of any kind. Wald/permutation-only.
   Houses: rank-based tests (Wilcoxon, sign test, Jonckheere-Terpstra, ridit), KM-based survival
   tests (log-rank, Gehan-Wilcox, RMST, KM-diff), exact tests, g-computation without a model
   likelihood, quantile regression (check-loss, not a true likelihood).

1. **`InferenceQuasiLik`** — an estimating-equation/working-variance construct; not a true density,
   so there is no valid likelihood-ratio statistic without ad hoc correction. Mechanically
   Wald-only, same as tier 0, but semantically distinct and a home for quasi-likelihood-specific
   shared logic (e.g. quasi-deviance helpers) as it gets built out.
   Houses: GEE (moving here fixes its current inconsistent wiring outside the whole subtree),
   quasi-Poisson, robust-Poisson, robust/sandwich-corrected regression.

2. **`InferencePartialLik`** — a genuine likelihood function with rigorously valid score/Wald/
   gradient/LR asymptotics (this is exactly what Cox 1975 proves), but it does not fully specify the
   joint density — a nuisance function is left unspecified (the baseline hazard for Cox, the
   conditioning distribution for conditional logistic).
   Houses: Cox PH, stratified Cox, conditional logistic regression, adjacent-category conditional
   logit, LWA-Cox.

3. **`InferenceFullLik`** — a fully specified generative density. Superset of tier 2's capability,
   plus trivial parametric simulation from `p(y|θ,X)`, plus eligibility for exact analytic
   corrections (Fisher information, Bartlett-exact) since the density is known in closed form.
   Houses: Poisson, NegBin, Hurdle/ZI mixtures, logistic, probit, cauchit, cloglog, ordered-probit,
   proportional-odds, stereotype-logit, log-binomial, beta, zero-one-inflated-beta, OLS, Weibull,
   Clayton-copula-with-parametric-margins, dependent-censoring transform.

### Mixins, each gated by a minimum tier, layered orthogonally

- **`InferenceMixinParamBootstrapSimulate`** — today's `InferenceParamBootstrap`, converted from an
  inheritance-chain class into an actual Pattern-1 mixin (like `InferenceMixinKKPassThrough` already
  is). Requires **≥ tier 2 (Partial)**. Provides `simulate_under_lik_null()` and the public
  `compute_lik_ratio_bootstrap_*` methods. A tier-2 class (Cox, via Breslow-plug-in semiparametric
  simulation) and a tier-3 class (Poisson, via direct parametric simulation) can both splice this in
  — which is exactly why Cox reaching "param bootstrap" was never actually a taxonomy violation: the
  mixin's requirement is "≥ Partial," not "== Full."
- **`InferenceMixinBartlettApprox`** — requires the simulate mixin above (already built this way
  today; unchanged).
- **`InferenceMixinBartlettExact`** / Cordeiro-Ferrari / Lemonte-gradient-approx — requires
  **== tier 3 (Full) only**. Gives the two currently-scaffolded-but-unused mixins
  (`InferenceMixinCordeiroFerrariApprox`, `InferenceMixinLemonteGradientApprox`) an actual home.
- **`InferenceMixinFisherInformation`** — requires **== tier 3 (Full) only**.
- **`InferenceMixinKKPassThroughCompound`** (already exists, unchanged) — orthogonal to all tiers.
  Caps the composing class's effective capability to Wald-only regardless of its components' tier,
  since a weighted combination of two separately-fit estimates has no single joint likelihood. Under
  this design every concrete IVWC compound class is `InferenceNoLik` + this mixin, regardless of
  whether its sub-model pieces are tier-2 (Cox) or tier-3 (Poisson) — composition, not likelihood
  strength, is what caps it.
- **`InferenceMixinKKPassThrough`** (already exists, already correctly a mixin) — unaffected.

### What happens to today's abstract backbone

| Today | Under this design |
|---|---|
| `InferenceAsymp`, `InferenceMLEorKMSummaryTable` | unchanged — pre-likelihood ancestors, no likelihood needed at any point above `InferenceAsympLik` |
| `InferenceAsympLik` | unchanged — stays the generic wald/score/gradient/lik_ratio/bartlett dispatcher; tiers 0–3 all sit under it |
| *(new)* `InferenceNoLik`, `InferenceQuasiLik`, `InferencePartialLik`, `InferenceFullLik` | replace the ~25+ scattered `supports_likelihood_tests = function() FALSE/TRUE` overrides with one authoritative definition per tier |
| `InferenceParamBootstrap` | retired as an inheritance-chain class; becomes `InferenceMixinParamBootstrapSimulate`, spliceable into tier 2 or tier 3 daughters |
| `InferenceCountLikelihood(NoParamBootstrap)` | plumbing (warm-starts, parameter packing) extracted into `InferenceMixinCountLikelihoodPlumbing`, spliced onto `InferenceFullLik` daughters; class retired |
| `InferenceAsympLikStdModCache(NoParamBootstrap)` | plumbing (caching) extracted into `InferenceMixinStdModCache`, spliced onto `InferenceFullLik` daughters; class retired — the single biggest mechanical change here, since this is the largest branch in the package (see the non-IVWC lattice) |
| `InferenceKKPassThroughCompound(NoParamBootstrap)` | retired as inheritance-chain classes; `InferenceMixinKKPassThroughCompound` gets spliced directly onto `InferenceNoLik` (every current concrete IVWC compound lands here, per the point above) |

### Where the original counterexamples land now

- **Cox reaching "param bootstrap"**: tier 2 + `InferenceMixinParamBootstrapSimulate` — the mixin's
  gate is "≥ Partial," so this was never a violation.
- **Quasi-Poisson flagged to zero capability**: correctly tier 1 (Quasi) — the taxonomy predicts
  exactly this, it isn't an anomaly.
- **GEE outside the `AsympLik` subtree**: a real inconsistency to fix by moving it into tier 1
  alongside quasi-Poisson — evidence the current *wiring* is incomplete, not that the taxonomy fails.
- **IVWC/CLMM compounds flagged to zero capability despite "full" sub-models**: `InferenceNoLik` +
  `InferenceMixinKKPassThroughCompound` — composition is orthogonal to tier, exactly your instinct
  that this should be a mixin.
- **`use_rcpp`-conditional capability drops** (9 classes across the two lattices): not a fifth tier —
  an engineering gap (the pure-R fallback was never built out to the same feature level as the Rcpp
  path). Flagged as follow-up tech debt, not modeled as a permanent axis.
- **Bartlett-exact implemented for exactly one class**: `InferenceMixinBartlettExact`, gated at
  tier 3 as necessary-but-not-sufficient — an optional top-tier mixin, not a fifth level.

The two published capability-lattice artifacts (`package_metadata/ivwc_capability_lattice.html`,
`package_metadata/non_ivwc_capability_lattice.html`) remain useful as the raw evidence base for this
mapping — every row in both tables now has an unambiguous (tier, mixin-set) home under this design.

### Scope and sequencing (decided)

**Decision: proceed with the six-file IVWC migration now, as already specified, and explicitly treat
those six files as a known future touch-point.** They land on `InferenceKKPassThroughCompoundNoParamBootstrap`
today; once this tier restructuring lands, that class is retired in favor of `InferenceNoLik` +
`InferenceMixinKKPassThroughCompound`, and the six files move a second time — a mechanical
`inherit =` + mixin-splice change, not a behavioral one, so the double-touch is cheap. This was chosen
over (a) blocking the six-file migration on the restructuring's implementation plan, which doesn't
exist yet, or (b) folding the six-file migration directly into the restructuring, which would leave a
ready, numerically-inert cleanup sitting idle for no benefit.

The restructuring itself remains out of scope for immediate implementation: converting
`InferenceParamBootstrap` into a mixin and retiring `InferenceAsympLikStdModCache` in particular touch
the majority of the package's ~130 kernels, so it needs its own dedicated implementation plan —
most plausibly executed tier-by-tier (introduce the four tier classes and the mixin conversions first,
with the existing abstract classes as deprecated pass-throughs, then migrate concrete families in
batches) rather than as one atomic change. That plan is not written yet.

## Follow-ups (explicitly not doing now)

- Auditing/removing the now-provably-dead `supports_lik_ratio_param_bootstrap` override on
  `InferenceAbstractKKWeibullFrailtyIVWC` (meaningless on the `NoParamBootstrap` branch).
- Any further consolidation of the three out-of-scope classes' smaller remaining duplication
  (their boilerplate doesn't fit the compound mold, so no shared base is proposed for them here).
