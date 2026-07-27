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

## Follow-ups (explicitly not doing now)

- Auditing/removing the now-provably-dead `supports_lik_ratio_param_bootstrap` override on
  `InferenceAbstractKKWeibullFrailtyIVWC` (meaningless on the `NoParamBootstrap` branch).
- Any further consolidation of the three out-of-scope classes' smaller remaining duplication
  (their boilerplate doesn't fit the compound mold, so no shared base is proposed for them here).
