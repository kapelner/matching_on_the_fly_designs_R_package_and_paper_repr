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
   FullLikelihood < ...`, so that likelihood taxonomy and capability level are the same axis encoded
   via inheritance depth. Rejected: the codebase already contains concrete counterexamples showing
   taxonomy and capability are orthogonal, not nested:
   - **Cox PH is partial-likelihood yet reaches top capability** (`StratCoxPHOneLik`, see #3 above) —
     a "partial" family outranking families the taxonomy would call more "full."
   - **Quasi-Poisson and robust-Poisson are quasi-likelihood yet reach zero capability** despite
     inheriting through `InferenceCountCompositeLikelihood -> InferenceParamBootstrap`: their
     `supports_likelihood_tests` and `supports_lik_ratio_param_bootstrap` are both hard-`FALSE` on
     the shared base class and neither subclass overrides them. Class-chain position and actual
     capability are decoupled.
   - **GEE, also quasi-likelihood, sits entirely outside the `InferenceAsympLik` subtree**
     (`inherit = InferenceAsymp` directly), with no likelihood-test dispatch machinery at all — while
     quasi-Poisson (same taxonomic bucket) is nested three levels inside that subtree. Two
     quasi-likelihood families are incomparable in capability, not merely unordered.
   - **`InferenceAbstractKKOrdinalCLMM` is a full joint (mixed-model) likelihood yet is flagged to
     zero capability** (`supports_likelihood_tests = FALSE`), same as the taxonomically "lesser"
     compound/IVWC classes.

   A linear chain requires every class to occupy exactly one rank in a total order; these families
   are genuinely incomparable on the taxonomy axis vs. the capability axis, which is why the
   codebase encodes capability as orthogonal boolean flags (`supports_likelihood_tests`,
   `supports_lik_ratio_param_bootstrap`, `supports_bartlett_likelihood_ratio_exact/approx`,
   `supports_fisher_information`) — several of them runtime-conditional functions rather than fixed
   per-class constants — dispatched generically from `InferenceAsympLik`, instead of via inheritance
   depth. That mechanism is left as-is.

   **Visual evidence — two capability-lattice matrices** (class × capability flag, saved as
   standalone self-contained HTML, viewable directly in a browser):
   `package_metadata/ivwc_capability_lattice.html` (the IVWC/OneLik pairs plus GEE, robust-regr,
   CLMM, paired-sign-test) and `package_metadata/non_ivwc_capability_lattice.html` (everything
   else: the 11-class abstract backbone plus every remaining concrete family across continuous,
   count, ordinal, proportion, incidence, and survival). The second sweep surfaced further
   confirmation beyond the four bullets above:

   - **A third major abstract branch**, `InferenceAsympLikStdModCache` / its `NoParamBootstrap`
     sibling — structurally the same role as `InferenceKKPassThroughCompound(NoParamBootstrap)` —
     turns out to be the single *largest* branch in the package (hosts most non-KK GLM-style
     families: logistic, probit, cauchit, cloglog, ordered-probit, proportional-odds,
     stereotype-logit, adjacent-category-logit, log-binomial, beta, zero-one-inflated-beta, Cox,
     stratified Cox, Weibull, dependent-censoring transform).
   - **The `use_rcpp`-conditional pattern recurs 8 more times** beyond Weibull-frailty: continuous/
     count/ordinal GLMM, non-KK Cox and stratified Cox, and the whole zero-augmented-Poisson family
     (Hurdle-Poisson/ZIP/ZINB) all silently drop to Wald-only on the pure-R fallback path — another
     case where capability is a runtime condition, not a chain position.
   - **A naming-convention break**: four classes named `...OneLik` (continuous and proportion
     quantile-regression, continuous robust-regression) are named like the full-capability branch
     established by Cox/Copula/Frailty/CondLogit/OLS, but actually inherit
     `InferenceKKPassThroughCompoundNoParamBootstrap` — zero likelihood-test capability.
   - **The "inside `InferenceParamBootstrap`'s static chain, flagged to zero capability" pattern**
     (first seen with quasi-Poisson/robust-Poisson) recurs on a much larger group: simple mean-diff
     (+ pooled-var), simple Wilcoxon, incidence Wald/CMH/extended-Robins, the KK-marginal g-comp
     family, fractional logit, risk-diff, and non-KK modified-Poisson — reinforcing that chain
     position and capability are decoupled throughout the package, not just in one or two spots.
   - A second **dead-flag instance** (`InferenceSurvivalKKWeibullMarginal`, outside the `AsympLik`
     subtree entirely, `supports_likelihood_tests = FALSE` but never consulted by anything) and a
     **terminology collision** (`InferenceOrdinalPartialProportionalOdds` — "partial" names the
     *partial proportional odds* statistical model, unrelated to "partial likelihood").

## Follow-ups (explicitly not doing now)

- Auditing/removing the now-provably-dead `supports_lik_ratio_param_bootstrap` override on
  `InferenceAbstractKKWeibullFrailtyIVWC` (meaningless on the `NoParamBootstrap` branch).
- Any further consolidation of the three out-of-scope classes' smaller remaining duplication
  (their boilerplate doesn't fit the compound mold, so no shared base is proposed for them here).
