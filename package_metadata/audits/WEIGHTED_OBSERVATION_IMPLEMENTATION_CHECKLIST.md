# Weighted-Observation Implementation Checklist

This file tracks which concrete bootstrap-capable inference paths currently have a
real implementation of
`compute_estimate_with_bootstrap_weights(subject_or_block_weights, ...)`
for the Bayesian bootstrap.

This audit is based on the current source tree, not on prior planning notes.
Inherited implementations from concrete or abstract family classes count as
implemented. The stub on `InferenceBayesianBootstrap` does not.

Per prior preference, the main tables below leave out `IVWC` paths. A short
`IVWC omitted` section is included separately at the end.

## Current Audit Summary

- In-scope concrete bootstrap-capable classes with real weighted support: `102`
  (`87` shown in the main table below, plus `13` IVWC classes and `2`
  Bai-adjusted classes that have a real *inherited surrogate* implementation
  but are kept out of the main table per prior preference — see "IVWC And
  Bai-Adjusted Classes Omitted From The Main Tables" below. This corrects a
  stale claim in the previous version of this file that these `15` classes
  were still unimplemented; see the re-audit note in that section.)
- In-scope concrete bootstrap-capable classes still unimplemented: `0`
- Exact/custom out-of-scope classes: `6`
- Non-IVWC unimplemented classes shown in the main tables below: `0`
- IVWC + Bai-adjusted classes omitted from the main tables (implemented via
  shared inherited surrogate, not because they lack an implementation): `15`

## Crosscheck Against The Previous Checklist

The previous version of this file was stale.

Main discrepancies:

- It said only `3` concrete paths were implemented. The current source has `84`.
- It marked these currently implemented classes as `Not implemented`:
  `InferenceAllSimpleMeanDiff`,
  `InferenceContinLin`,
  `InferenceContinOLS`,
  `InferenceCountPoisson`,
  `InferenceCountPoissonKKGEE`,
  `InferenceIncidKKGEE`,
  `InferenceIncidLogRegr`,
  `InferenceIncidModifiedPoisson`,
  `InferenceIncidRiskDiff`,
  `InferenceOrdinalGCompMeanDiff`,
  `InferenceOrdinalKKGEE`,
  `InferenceOrdinalPropOddsRegr`,
  `InferencePropKKGEE`,
  `InferenceSurvivalGehanWilcox`,
  `InferenceSurvivalKMDiff`,
  `InferenceSurvivalLogRank`,
  `InferenceSurvivalRestrictedMeanDiff`.
- It omitted these concrete classes entirely:
  `InferenceOrdinalStereotypeLogitRegr`,
  `InferenceOrdinalContRatioRegr`,
.
- It listed custom extension classes in the tracked table even though they should
  remain out of scope for the package’s built-in weighted rollout.

## Implemented

These classes currently have a real weighted implementation, either directly or
through a non-stub parent class.

Two additional implemented paths are omitted from the main table because they
are `IVWC`:

- `InferenceCountKKHurdlePoissonIVWC`

### Implemented Non-IVWC Paths

| Class | Family | Provider | Notes |
|---|---|---|---|
| `InferenceAllSimpleMeanDiff` | All / continuous | `InferenceAllSimpleMeanDiff` | Weighted empirical mean difference |
| `InferenceAllSimpleMeanDiffPooledVar` | All / continuous | `InferenceAllSimpleMeanDiff` | Inherits weighted mean-difference estimate |
| `InferenceAllSimpleWilcox` | All / nonparametric | `InferenceAllSimpleWilcox` | Weighted Hodges-Lehmann surrogate via weighted pairwise-difference median |
| `InferenceContinLin` | Continuous | `InferenceContinLin` | Weighted `lm.wfit` on Lin design |
| `InferenceContinKKGLMM` | Continuous KK GLMM | `InferenceMixinKKGLMMShared` | Weighted `glmmTMB` random-intercept surrogate |
| `InferenceContinKKOLSOneLik` | Continuous KK one-likelihood | `InferenceContinKKOLSOneLik` | Weighted stacked `lm.wfit` |
| `InferenceContinKKQuantileRegrOneLik` | Continuous KK one-likelihood | `InferenceAbstractKKQuantileRegrOneLik` | Weighted stacked `quantreg::rq` |
| `InferenceContinKKRobustRegrOneLik` | Continuous KK one-likelihood | `InferenceAbstractKKRobustRegrOneLik` | Weighted stacked robust-regression surrogate |
| `InferenceContinOLS` | Continuous | `InferenceContinOLS` | Weighted `lm.wfit` |
| `InferenceContinQuantileRegr` | Continuous | `InferenceContinQuantileRegr` | Weighted `quantreg::rq` |
| `InferenceContinRobustRegr` | Continuous | `InferenceContinRobustRegr` | Weighted `MASS::rlm` |
| `InferenceCountCompositeLikelihood` | Count composite likelihood | `InferenceCountCompositeLikelihood` | Weighted companion-Poisson surrogate |
| `InferenceCountKKCondPoissonOneLik` | Count KK one-likelihood | `InferenceAbstractKKPoissonCondPoissonOneLik` | Weighted one-likelihood path |
| `InferenceCountKKGLMM` | Count KK GLMM | `InferenceMixinKKGLMMShared` | Weighted `glmmTMB` Poisson mixed-model surrogate |
| `InferenceCountKKHurdlePoissonOneLik` | Count KK one-likelihood | `InferenceAbstractKKHurdlePoissonOneLik` | Weighted one-likelihood path |
| `InferenceCountNegBin` | Count | `InferenceCountNegBin` | Weighted `MASS::glm.nb` |
| `InferenceCountHurdleNegBin` | Count | `InferenceCountHurdleNegBin` | Weighted `glmmTMB` hurdle negative binomial |
| `InferenceCountPoisson` | Count | `InferenceCountPoisson` | First-wave weighted Poisson |
| `InferenceCountPoissonKKGEE` | Count KK GEE | `InferenceCountPoissonKKGEE` | Native weighted KK GEE |
| `InferenceCountQuasiPoisson` | Count | `InferenceCountQuasiPoisson` | Weighted Poisson point refit for quasi-Poisson estimate |
| `InferenceCountRobustPoisson` | Count | `InferenceCountRobustPoisson` | Weighted Poisson point refit for robust Poisson estimate |
| `InferenceCountHurdlePoisson` | Count | `InferenceCountZeroAugmentedPoissonAbstract` | Weighted `glmmTMB` hurdle Poisson |
| `InferenceCountZeroInflatedPoisson` | Count | `InferenceCountZeroAugmentedPoissonAbstract` | Weighted `glmmTMB` zero-inflated Poisson |
| `InferenceCountZeroInflatedNegBin` | Count | `InferenceCountZeroAugmentedPoissonAbstract` | Weighted `glmmTMB` zero-inflated negative binomial |
| `InferenceIncidCMH` | Incidence | `InferenceAllSimpleMeanDiff` | Weighted empirical contrast path |
| `InferenceIncidBinomialIdentityRiskDiff` | Incidence | `InferenceIncidBinomialIdentityRiskDiff` | Weighted constrained identity-binomial |
| `InferenceIncidExtendedRobins` | Incidence | `InferenceAllSimpleMeanDiff` | Weighted empirical contrast path |
| `InferenceIncidGCompRiskDiff` | Incidence g-comp | `InferenceIncidGCompAbstract` | Weighted g-computation |
| `InferenceIncidGCompRiskRatio` | Incidence g-comp | `InferenceIncidGCompAbstract` | Weighted g-computation |
| `InferenceIncidKKGEE` | Incidence KK GEE | `InferenceIncidKKGEE` | Native weighted KK GEE |
| `InferenceIncidKKCondLogitOneLik` | Incidence KK one-likelihood | `InferenceAbstractKKCondLogitOneLik` | Weighted combined logistic surrogate |
| `InferenceIncidKKCondLogitGLMMOneLik` | Incidence KK hybrid one-likelihood | `InferenceAbstractKKCondLogitGLMM` | Weighted logistic mixed-model surrogate |
| `InferenceIncidKKGCompRiskDiff` | Incidence KK g-comp | `InferenceIncidKKGCompAbstract` | Weighted KK g-computation |
| `InferenceIncidKKGCompRiskRatio` | Incidence KK g-comp | `InferenceIncidKKGCompAbstract` | Weighted KK g-computation |
| `InferenceIncidKKModifiedPoisson` | Incidence KK marginal | `InferenceAbstractKKModifiedPoisson` | Weighted KK modified-Poisson point refit |
| `InferenceIncidKKNewcombeRiskDiff` | Incidence KK non-likelihood | `InferenceIncidKKNewcombeRiskDiff` | Weighted empirical risk-difference surrogate |
| `InferenceIncidLogBinomial` | Incidence | `InferenceIncidLogBinomial` | Weighted constrained log-binomial |
| `InferenceIncidLogRegr` | Incidence | `InferenceIncidLogRegr` | First-wave weighted logistic |
| `InferenceIncidModifiedPoisson` | Incidence | `InferenceIncidModifiedPoisson` | First-wave weighted modified Poisson |
| `InferenceIncidMiettinenNurminenRiskDiff` | Incidence | `InferenceIncidMiettinenNurminenRiskDiff` | Weighted empirical surrogate for MN bootstrap estimate |
| `InferenceIncidNewcombeRiskDiff` | Incidence | `InferenceIncidNewcombeRiskDiff` | Weighted empirical Newcombe point estimate |
| `InferenceIncidProbitRegr` | Incidence | `InferenceIncidProbitRegr` | Weighted `glm.fit` probit surrogate |
| `InferenceIncidRiskDiff` | Incidence | `InferenceIncidRiskDiff` | Weighted identity-scale risk difference |
| `InferenceIncidWald` | Incidence | `InferenceAllSimpleMeanDiff` | Weighted empirical contrast path |
| `InferenceOrdinalGCompMeanDiff` | Ordinal g-comp | `InferenceOrdinalGCompMeanDiff` | Weighted ordinal g-computation |
| `InferenceOrdinalAdjCatLogitRegr` | Ordinal | `InferenceOrdinalAdjCatLogitRegr` | Weighted cumulative-logit surrogate for adjacent-category estimate |
| `InferenceOrdinalCloglogRegr` | Ordinal | `InferenceOrdinalCloglogRegr` | Weighted cloglog ordinal surrogate |
| `InferenceOrdinalCauchitRegr` | Ordinal | `InferenceOrdinalCauchitRegr` | Weighted cauchit ordinal surrogate |
| `InferenceOrdinalContRatioRegr` | Ordinal | `InferenceOrdinalContRatioRegr` | Weighted ordinal surrogate for continuation-ratio estimate |
| `InferenceOrdinalOrderedProbitRegr` | Ordinal | `InferenceOrdinalOrderedProbitRegr` | Weighted ordered-probit surrogate |
| `InferenceOrdinalKKGEE` | Ordinal KK GEE | `InferenceOrdinalKKGEE` | Weighted ordinal surrogate path |
| `InferenceOrdinalKKCLMM` | Ordinal KK CLMM | `InferenceAbstractKKOrdinalCLMM` | Weighted cumulative-link mixed-model surrogate with native weighted proportional-odds fallback and ordinal-link surrogates |
| `InferenceOrdinalKKCLMMProbit` | Ordinal KK CLMM | `InferenceAbstractKKOrdinalCLMM` | Weighted cumulative-link mixed-model surrogate with probit-link ordinal surrogate |
| `InferenceOrdinalKKCLMMCauchit` | Ordinal KK CLMM | `InferenceAbstractKKOrdinalCLMM` | Weighted cumulative-link mixed-model surrogate with cauchit-link ordinal surrogate |
| `InferenceOrdinalKKCLMMCloglog` | Ordinal KK CLMM | `InferenceAbstractKKOrdinalCLMM` | Weighted cumulative-link mixed-model surrogate with cloglog-link ordinal surrogate |
| `InferenceOrdinalKKCondAdjCatLogitRegr` | Ordinal KK likelihood | `InferenceOrdinalKKCondAdjCatLogitRegr` | Weighted ordinal surrogate |
| `InferenceOrdinalKKGLMM` | Ordinal KK GLMM | `InferenceMixinKKGLMMShared` | Weighted cumulative-link mixed-model surrogate |
| `InferenceOrdinalJonckheereTerpstraTest` | Ordinal | `InferenceOrdinalJonckheereTerpstraTest` | Weighted empirical superiority / JT surrogate |
| `InferenceOrdinalPartialProportionalOddsRegr` | Ordinal | `InferenceOrdinalPartialProportionalOddsRegr` | Weighted PPO with native weighted proportional-odds fast path and weighted VGAM/CLM/POLR fallbacks |
| `InferenceOrdinalRidit` | Ordinal | `InferenceOrdinalRidit` | Weighted empirical ridit analysis |
| `InferenceOrdinalPropOddsRegr` | Ordinal | `InferenceOrdinalPropOddsRegr` | Weighted proportional odds |
| `InferenceOrdinalStereotypeLogitRegr` | Ordinal | `InferenceOrdinalStereotypeLogitRegr` | Weighted cumulative-logit surrogate for stereotype-logit estimate |
| `InferenceOrdinalPairedSignTest` | Ordinal | `InferenceOrdinalPairedSignTest` | Weighted matched-pair sign statistic |
| `InferencePropBetaRegr` | Proportion | `InferencePropBetaRegr` | Weighted beta-regression / weighted logit fallback |
| `InferencePropFractionalLogit` | Proportion | `InferencePropFractionalLogit` | Weighted fractional logit |
| `InferencePropZeroOneInflatedBetaRegr` | Proportion | `InferencePropZeroOneInflatedBetaRegr` | Weighted mixture surrogate with weighted beta/logit components |
| `InferencePropGCompMeanDiff` | Proportion g-comp | `InferencePropGCompAbstract` | Weighted g-computation |
| `InferencePropKKGEE` | Proportion KK GEE | `InferencePropKKGEE` | Native weighted KK GEE |
| `InferencePropKKGLMM` | Proportion KK GLMM | `InferenceAbstractKKLogisticGLMMOneLik` | Weighted logistic mixed-model surrogate |
| `InferencePropKKQuantileRegrOneLik` | Proportion KK one-likelihood | `InferenceAbstractKKQuantileRegrOneLik` | Weighted stacked `quantreg::rq` on transformed response |
| `InferenceSurvivalGehanWilcox` | Survival | `InferenceSurvivalGehanWilcox` | Weighted Peto-Prentice score estimate |
| `InferenceSurvivalCoxPHRegr` | Survival | `InferenceSurvivalCoxPHRegr` | Weighted Cox PH surrogate via `survival::coxph` |
| `InferenceSurvivalDepCensTransformRegr` | Survival | `InferenceSurvivalDepCensTransformRegr` | Weighted Cox-style surrogate for dependent-censoring transform estimate |
| `InferenceSurvivalKMDiff` | Survival | `InferenceSurvivalKMDiff` | Weighted KM median difference |
| `InferenceSurvivalKKClaytonCopulaOneLik` | Survival KK one-likelihood | `InferenceAbstractKKClaytonCopulaOneLik` | Weighted Weibull AFT surrogate over combined KK data |
| `InferenceSurvivalKKLWACoxPHOneLik` | Survival KK one-likelihood | `InferenceAbstractKKLWACoxOneLik` | Weighted marginal Cox PH surrogate with KK clusters |
| `InferenceSurvivalKKStratCoxPHOneLik` | Survival KK one-likelihood | `InferenceAbstractKKStratCoxOneLik` | Weighted stratified Cox PH surrogate over KK units |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | Survival KK one-likelihood | `InferenceAbstractKKWeibullFrailtyOneLik` | Weighted Weibull AFT surrogate with KK clusters |
| `InferenceSurvivalLogRank` | Survival | `InferenceSurvivalLogRank` | Weighted log-rank score estimate |
| `InferenceSurvivalRestrictedMeanDiff` | Survival | `InferenceSurvivalRestrictedMeanDiff` | Weighted RMST difference |
| `InferenceSurvivalStratCoxPHRegr` | Survival | `InferenceSurvivalStratCoxPHRegr` | Weighted stratified-Cox surrogate |
| `InferenceSurvivalWeibullRegr` | Survival | `InferenceSurvivalWeibullRegr` | Weighted Weibull AFT surrogate via `survival::survreg` |

## Still Unimplemented

There are no remaining non-IVWC concrete rollout targets that still fall
through to the Bayesian-bootstrap stub, and — per the re-audit below — there
are no remaining unimplemented concrete rollout targets at all. The IVWC
block and the Bai-adjusted KK compound family (`InferenceBaiAdjustedTKK14`,
`InferenceBaiAdjustedTKK21`) are omitted from the main tables per prior
preference, but both actually inherit a real (if surrogate/approximate)
`compute_estimate_with_bootstrap_weights` implementation — see the corrected
section below. The previous version of this checklist claimed the
Bai-adjusted family was "explicitly gated off from Bayesian bootstrap"; that
is no longer accurate (and may never have been enforced in code — see
TODO-1).

### Near-Term / Moderate Difficulty

These are the most plausible next implementations if we continue expanding
coverage without adding major new infrastructure.

| Class | Family | Why It Is Not Done Yet |
|---|---|---|
No remaining classes in this moderate-difficulty block after the current implementation tranche.

### Ordinal Likelihood Gap Block

The previous ordinal-likelihood gap block is now cleared for the current
rollout. The implementations here are a mix of exact weighted empirical
definitions (`Ridit`, `PairedSignTest`) and documented weighted surrogate fits
for ordinal link families whose native weighted backends do not yet exist.

### KK One-Likelihood / GLMM / Compound Paths

The current rollout now covers this block with a mix of weighted stacked-design
re-estimation, weighted `glmmTMB` mixed-model surrogates, and direct empirical
surrogates where the original KK estimating equation did not yet have a native
weighted backend.

### Survival Regression / Survival KK Deferred Paths

The previous deferred survival block is now cleared for the current rollout.
These implementations use weighted survival-model surrogates rather than native
weighted versions of the original estimating procedures.

## IVWC And Bai-Adjusted Classes Omitted From The Main Tables

**Re-audit correction (this pass):** the previous version of this file
described these classes as "still unimplemented." That is stale. All of them
inherit `InferenceKKPassThroughCompoundNoParamBootstrap`
(`EDI/R/inference_all_abstract_KK_passthrough_compound.R:108`), whose public
component (`inference_kk_passthrough_compound_components`, composed at
`EDI/R/inference_all_abstract_KK_passthrough_compound.R:8-89`) supplies a real
`compute_estimate_with_bootstrap_weights` (line 36): it computes a
matched-pair-difference / reservoir-mean-difference surrogate combined via the
observed-fit `w_star` weight (lines 48-86). This is the same "surrogate"
pattern already credited as `Implemented` for ~30 other rows in the main
table above (e.g. `InferenceIncidCMH`, `InferenceIncidWald`, all via
`InferenceAllSimpleMeanDiff`'s "Weighted empirical contrast path"). By this
document's own stated rule ("Inherited implementations from concrete or
abstract family classes count as implemented"), these classes are
implemented. They remain **excluded from the main table above purely for
organizational/display preference** (to keep IVWC paths out of the main
inventory), not because they lack an implementation.

The class-specific override hook this fallback checks for,
`private$compute_weighted_estimate_ivwc` (line 44), is currently defined by
**zero** concrete classes (`grep -rln "compute_weighted_estimate_ivwc\s*=\s*function" EDI/R/*.R`
returns nothing) — every class below runs the generic surrogate fallback, not
a class-specific weighted re-fit. See TODO-2.

- `InferenceAllKKMeanDiffIVWC`
- `InferenceAllKKWilcoxIVWC`
- `InferenceContinKKOLSIVWC`
- `InferenceContinKKQuantileRegrIVWC`
- `InferenceContinKKRobustRegrIVWC`
- `InferenceIncidKKCondLogitIVWC`
- `InferenceIncidKKCondLogitGLMMIVWC`
- `InferencePropKKQuantileRegrIVWC`
- `InferenceSurvivalKKClaytonCopulaIVWC`
- `InferenceSurvivalKKLWACoxPHIVWC`
- `InferenceSurvivalKKRankRegrIVWC`
- `InferenceSurvivalKKStratCoxPHIVWC`
- `InferenceSurvivalKKWeibullFrailtyIVWC`
- `InferenceBaiAdjustedTKK14` (`EDI/R/inference_continuous_KK14_bai.R:19`,
  via `InferenceBaiAdjustedT` → `InferenceKKPassThroughCompoundNoParamBootstrap`)
- `InferenceBaiAdjustedTKK21` (`EDI/R/inference_continuous_KK21_bai.R:19`,
  same inheritance path) — neither Bai class was listed in any table
  (`Implemented`, `Still Unimplemented`, `IVWC Omitted`, or `Out Of Scope`) in
  the previous version of this file; this was a bookkeeping gap, now closed.

## Out Of Scope

These are not part of the package’s built-in weighted Bayesian-bootstrap rollout.

- `InferenceIncidExactFisher`
- `InferenceIncidExactBinomial`
- `InferenceIncidenceExactZhang`
- `InferenceCustomAsymp`
- `InferenceCustomRand`
- `InferenceCustomBoot`

## Recommended Implementation Order For The Near-Term Block

This is the recommended implementation order for the still-unimplemented
near-term families, based on backend availability, expected semantic clarity
under Bayesian-bootstrap weights, and numerical risk.

The previous near-term, ordinal-gap, KK one-likelihood / GLMM compound, and
survival-regression blocks have been cleared. The remaining in-scope work is
limited to any future refinement of native weighted backends for the IVWC and
Bai-adjusted paths (see "IVWC And Bai-Adjusted Classes Omitted From The Main
Tables" above) — those paths already have a real surrogate implementation, so
this is a refinement, not a rollout gap.

## Suggested Milestones

The old near-term and ordinal milestones are complete.

## Backend Selection Rule

- Prefer native backend extension where the class already uses a package-owned
  optimizer.
- Use external weighted fits only when the class already depends on that
  external engine or when the implementation is clearly temporary and
  documented as such.
- Keep `estimate_only = TRUE` as the first target for every class, and defer
  weighted replicate-SE logic unless a Bayesian-bootstrap summary actually
  needs it.

## Implementation TODOs

### Dead Gate Functions (Correctness / Cleanup)

- [ ] TODO-1: `stop_bayesian_bootstrap_for_bai()` and
  `stop_bayesian_bootstrap_for_ivwc()` (`EDI/R/globals.R:89-105`) are defined
  but never called anywhere in `EDI/R/*.R` (confirmed via
  `grep -rn "stop_bayesian_bootstrap_for_bai\|stop_bayesian_bootstrap_for_ivwc" EDI/R/*.R`
  — zero call sites). They read as leftover guards from a version of the code
  where Bai/IVWC paths genuinely raised on Bayesian-bootstrap weighted
  re-estimation; today those paths silently run the generic surrogate
  fallback instead (`EDI/R/inference_all_abstract_KK_passthrough_compound.R:36-87`).
  Decide and act: either (a) delete both now-dead functions, or (b) if the
  surrogate fallback was never meant to be a real substitute for Bai/IVWC
  paths, wire them back in and treat those `15` classes as genuinely
  unimplemented again. Current code behavior (silent surrogate, no warning)
  suggests (a), but this is a product decision, not just a doc fix.

### Class-Specific Weighted Refits (Backend Refinement)

- [ ] TODO-2: Zero concrete classes currently define
  `private$compute_weighted_estimate_ivwc` (confirmed via
  `grep -rln "compute_weighted_estimate_ivwc\s*=\s*function" EDI/R/*.R` —
  empty result), so every one of the 15 classes listed in "IVWC And
  Bai-Adjusted Classes Omitted From The Main Tables" runs the same generic
  matched-difference/reservoir-mean surrogate
  (`EDI/R/inference_all_abstract_KK_passthrough_compound.R:48-86`) regardless
  of its actual estimator (OLS, robust regression, quantile regression, Cox
  PH, Clayton copula, Weibull frailty, Bai-adjusted T, ...). Per the
  "Backend Selection Rule" above, add class-specific
  `compute_weighted_estimate_ivwc` overrides for the highest-traffic families
  first — `InferenceContinKKOLSIVWC` and `InferenceContinKKRobustRegrIVWC`
  are the most natural first targets since their non-IVWC one-likelihood
  siblings (`InferenceContinKKOLSOneLik`, `InferenceContinKKRobustRegrOneLik`)
  already have dedicated weighted implementations per the main table above —
  use those as the reference implementation to adapt.
- [ ] TODO-3: The Bai-adjusted classes (`InferenceBaiAdjustedTKK14`,
  `InferenceBaiAdjustedTKK21`) use the same generic surrogate despite Bai
  adjustment being a distinct statistical correction from the plain
  matched-pair/reservoir estimate it falls back to. Verify whether the
  surrogate is a statistically valid stand-in for the Bai-adjusted estimator
  under Bayesian-bootstrap weights, or whether a Bai-specific weighted
  estimate is needed — this is a correctness question, not just a coverage
  gap, since TODO-1 already establishes nothing currently prevents this path
  from running and returning a result that looks like a real Bai estimate.

### Documentation Maintenance

- [ ] TODO-4: Keep this file's "Current Audit Summary" counts and the
  `Implemented` / `IVWC And Bai-Adjusted Classes Omitted` sections in sync
  going forward — re-run the two greps in TODO-1/TODO-2 whenever a class in
  the KK-passthrough-compound family changes, since this pass found the
  previous version of this file had silently drifted (claiming 15
  implemented classes were "still unimplemented"/"gated off") without any
  corresponding code change being the cause — the doc simply hadn't been
  re-verified against the code it was auditing.
