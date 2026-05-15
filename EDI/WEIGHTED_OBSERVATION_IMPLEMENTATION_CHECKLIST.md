# Weighted-Observation Implementation Checklist

This file tracks which concrete inference paths currently implement
`compute_estimate_with_bootstrap_weights(subject_or_block_weights, ...)`
for the Bayesian bootstrap.

Status meaning:

- `Implemented`: the class defines weighted re-estimation explicitly.
- `Not implemented`: the class does not currently define the weighted hook.
- `Out of scope`: exact-only, randomization-only, or custom extension paths that
  are not part of the current weighted Bayesian-bootstrap rollout.

Implementation status below is based mechanically on whether the class defines
`compute_estimate_with_bootstrap_weights(...)` in the current codebase.

## Implemented

| Class | Family | Path type | Status | Notes |
|---|---|---|---|---|
| `InferenceIncidLogRegr` | Incidence | Likelihood GLM | Implemented | First-wave weighted logistic support |
| `InferenceCountPoisson` | Count | Likelihood GLM | Implemented | First-wave weighted Poisson support |
| `InferenceIncidModifiedPoisson` | Incidence | Likelihood GLM | Implemented | First-wave weighted modified Poisson support |

## Continuous and All-Outcome Non-Likelihood Paths

| Class | Family | Path type | Status | Notes |
|---|---|---|---|---|
| `InferenceAllSimpleMeanDiff` | All | Non-likelihood | Not implemented | Simple mean-difference path |
| `InferenceAllSimpleMeanDiffPooledVar` | All | Non-likelihood | Not implemented | Pooled-variance mean-difference path |
| `InferenceAllSimpleWilcox` | All | Non-likelihood | Not implemented | Wilcoxon path |
| `InferenceContinLin` | Continuous | Non-likelihood | Not implemented | Linear-model contrast path |
| `InferenceContinOLS` | Continuous | Non-likelihood | Not implemented | OLS path |
| `InferenceContinRobustRegr` | Continuous | Non-likelihood | Not implemented | Robust regression path |
| `InferenceContinQuantileRegr` | Continuous | Non-likelihood | Not implemented | Quantile-regression path |
| `InferenceBaiAdjustedT` | Continuous | Non-likelihood | Not implemented | Bai-adjusted t path |
| `InferenceBaiAdjustedTKK14` | Continuous | Non-likelihood | Not implemented | KK14 Bai-adjusted t path |
| `InferenceBaiAdjustedTKK21` | Continuous | Non-likelihood | Not implemented | KK21 Bai-adjusted t path |
| `InferenceAllKKMeanDiffIVWC` | Continuous | KK non-likelihood | Not implemented | KK IV-weighted mean difference |
| `InferenceAllKKWilcoxIVWC` | Continuous | KK non-likelihood | Not implemented | KK IV-weighted Wilcoxon |
| `InferenceContinKKOLSIVWC` | Continuous | KK non-likelihood / IVWC | Not implemented | KK OLS IVWC path |
| `InferenceContinKKOLSOneLik` | Continuous | KK one-likelihood | Not implemented | KK one-likelihood OLS path |
| `InferenceContinKKRobustRegrIVWC` | Continuous | KK non-likelihood / IVWC | Not implemented | KK robust regression IVWC |
| `InferenceContinKKRobustRegrOneLik` | Continuous | KK one-likelihood | Not implemented | KK robust regression one-likelihood |
| `InferenceContinKKQuantileRegrIVWC` | Continuous | KK non-likelihood / IVWC | Not implemented | KK quantile regression IVWC |
| `InferenceContinKKQuantileRegrOneLik` | Continuous | KK one-likelihood | Not implemented | KK quantile regression one-likelihood |
| `InferenceContinKKGLMM` | Continuous | KK GLMM | Not implemented | Mixed-model path |

## Incidence Paths

| Class | Family | Path type | Status | Notes |
|---|---|---|---|---|
| `InferenceIncidenceWald` | Incidence | Non-likelihood | Not implemented | Wald path |
| `InferenceIncidCMH` | Incidence | Non-likelihood | Not implemented | CMH path |
| `InferenceIncidExtendedRobins` | Incidence | Non-likelihood | Not implemented | Extended Robins path |
| `InferenceIncidNewcombeRiskDiff` | Incidence | Non-likelihood | Not implemented | Newcombe risk difference |
| `InferenceIncidMiettinenNurminenRiskDiff` | Incidence | Non-likelihood | Not implemented | Miettinen-Nurminen risk difference |
| `InferenceIncidRiskDiff` | Incidence | Likelihood-style | Not implemented | Binomial risk-difference model path |
| `InferenceIncidBinomialIdentityRiskDiff` | Incidence | Likelihood GLM | Not implemented | Binomial identity-link path |
| `InferenceIncidLogBinomial` | Incidence | Likelihood GLM | Not implemented | Log-binomial path |
| `InferenceIncidProbitRegr` | Incidence | Likelihood GLM | Not implemented | Probit path |
| `InferenceIncidGCompRiskDiff` | Incidence | G-computation | Not implemented | Non-likelihood g-comp path |
| `InferenceIncidGCompRiskRatio` | Incidence | G-computation | Not implemented | Non-likelihood g-comp path |
| `InferenceIncidKKNewcombeRiskDiff` | Incidence | KK non-likelihood / IVWC | Not implemented | KK Newcombe path |
| `InferenceIncidKKGCompRiskDiff` | Incidence | KK g-computation | Not implemented | Appears in KK marginal file |
| `InferenceIncidKKGCompRiskRatio` | Incidence | KK g-computation | Not implemented | Appears in KK marginal file |
| `InferenceIncidKKModifiedPoisson` | Incidence | KK marginal / likelihood-style | Not implemented | KK modified Poisson |
| `InferenceIncidKKClogitIVWC` | Incidence | KK IVWC likelihood | Not implemented | Conditional logistic IVWC |
| `InferenceIncidKKClogitOneLik` | Incidence | KK one-likelihood | Not implemented | Conditional logistic one-likelihood |
| `InferenceIncidKKClogitPlusGLMMIVWC` | Incidence | KK IVWC hybrid | Not implemented | Clogit + GLMM IVWC |
| `InferenceIncidKKClogitPlusGLMMOneLik` | Incidence | KK one-likelihood hybrid | Not implemented | Clogit + GLMM one-likelihood |
| `InferenceIncidKKGEE` | Incidence | KK GEE | Not implemented | Explicit GEE path |
| `InferenceIncidKKGLMM` | Incidence | KK GLMM | Not implemented | Mixed-model path |

## Count Paths

| Class | Family | Path type | Status | Notes |
|---|---|---|---|---|
| `InferenceCountNegBin` | Count | Likelihood | Not implemented | Negative binomial |
| `InferenceCountQuasiPoisson` | Count | Likelihood-style | Not implemented | Quasi-Poisson |
| `InferenceCountRobustPoisson` | Count | Composite / robust likelihood | Not implemented | Sandwich-based Poisson |
| `InferenceCountCompositeLikelihood` | Count | Composite likelihood | Not implemented | Abstracted count composite path |
| `InferenceCountHurdlePoisson` | Count | Likelihood | Not implemented | Hurdle Poisson |
| `InferenceCountHurdleNegBin` | Count | Likelihood | Not implemented | Hurdle negative binomial |
| `InferenceCountZeroInflatedPoisson` | Count | Likelihood | Not implemented | Zero-inflated Poisson |
| `InferenceCountZeroInflatedNegBin` | Count | Likelihood | Not implemented | Zero-inflated negative binomial |
| `InferenceCountKKCPoissonIVWC` | Count | KK IVWC likelihood | Not implemented | Compound Poisson IVWC |
| `InferenceCountKKCPoissonOneLik` | Count | KK one-likelihood | Not implemented | Compound Poisson one-likelihood |
| `InferenceCountKKHurdlePoissonIVWC` | Count | KK IVWC likelihood | Not implemented | KK hurdle Poisson IVWC |
| `InferenceCountKKHurdlePoissonOneLik` | Count | KK one-likelihood | Not implemented | KK hurdle Poisson one-likelihood |
| `InferenceCountKKGLMM` | Count | KK GLMM | Not implemented | Mixed-model path |
| `InferenceCountPoissonKKGEE` | Count | KK GEE | Not implemented | Explicit GEE path |
| `InferenceCountPoissonUnivKKGEE` | Count | KK GEE | Not implemented | Univariable KK GEE |
| `InferenceCountPoissonMultiKKGEE` | Count | KK GEE | Not implemented | Multivariable KK GEE |

## Proportion Paths

| Class | Family | Path type | Status | Notes |
|---|---|---|---|---|
| `InferencePropFractionalLogit` | Proportion | Likelihood | Not implemented | Fractional logit |
| `InferencePropBetaRegr` | Proportion | Likelihood | Not implemented | Beta regression |
| `InferencePropZeroOneInflatedBetaRegr` | Proportion | Likelihood | Not implemented | Zero/one-inflated beta |
| `InferencePropGCompMeanDiff` | Proportion | G-computation | Not implemented | Non-likelihood g-comp path |
| `InferencePropKKGEE` | Proportion | KK GEE | Not implemented | Explicit GEE path |
| `InferencePropKKGLMM` | Proportion | KK GLMM | Not implemented | Mixed-model path |
| `InferencePropKKQuantileRegrIVWC` | Proportion | KK IVWC | Not implemented | Quantile-regression-based KK path |
| `InferencePropKKQuantileRegrOneLik` | Proportion | KK one-likelihood | Not implemented | Quantile-regression-based KK path |

## Ordinal Paths

| Class | Family | Path type | Status | Notes |
|---|---|---|---|---|
| `InferenceOrdinalRidit` | Ordinal | Non-likelihood | Not implemented | Ridit path |
| `InferenceOrdinalPairedSignTest` | Ordinal | Likelihood-style | Not implemented | Paired sign-test path |
| `InferenceOrdinalPropOddsRegr` | Ordinal | Likelihood | Not implemented | Proportional-odds |
| `InferenceOrdinalAdjCatLogitRegr` | Ordinal | Likelihood | Not implemented | Adjacent-category logit |
| `InferenceOrdinalCloglogRegr` | Ordinal | Likelihood | Not implemented | Cloglog ordinal |
| `InferenceOrdinalCauchitRegr` | Ordinal | Likelihood | Not implemented | Cauchit ordinal |
| `InferenceOrdinalOrderedProbitRegr` | Ordinal | Likelihood | Not implemented | Ordered probit |
| `InferenceOrdinalGCompMeanDiff` | Ordinal | G-computation | Not implemented | Non-likelihood g-comp path |
| `InferenceOrdinalKKGEE` | Ordinal | KK GEE | Not implemented | Explicit GEE path |
| `InferenceOrdinalKKGLMM` | Ordinal | KK GLMM | Not implemented | Mixed-model path |
| `InferenceOrdinalKKCondAdjCatLogitRegr` | Ordinal | KK likelihood | Not implemented | Conditional adjacent-category logit |
| `InferenceOrdinalKKCLMM` | Ordinal | KK CLMM | Not implemented | Cumulative link mixed model |
| `InferenceOrdinalKKCLMMProbit` | Ordinal | KK CLMM | Not implemented | Probit CLMM |
| `InferenceOrdinalKKCLMMCauchit` | Ordinal | KK CLMM | Not implemented | Cauchit CLMM |
| `InferenceOrdinalKKCLMMCloglog` | Ordinal | KK CLMM | Not implemented | Cloglog CLMM |

## Survival Paths

| Class | Family | Path type | Status | Notes |
|---|---|---|---|---|
| `InferenceSurvivalKMDiff` | Survival | Non-likelihood | Not implemented | Kaplan-Meier difference |
| `InferenceSurvivalRestrictedMeanDiff` | Survival | Non-likelihood | Not implemented | RMST path |
| `InferenceSurvivalLogRank` | Survival | Non-likelihood | Not implemented | Log-rank |
| `InferenceSurvivalGehanWilcox` | Survival | Non-likelihood | Not implemented | Gehan-Wilcoxon |
| `InferenceSurvivalCoxPHRegr` | Survival | Likelihood | Not implemented | Cox PH |
| `InferenceSurvivalStratCoxPHRegr` | Survival | Likelihood | Not implemented | Stratified Cox PH |
| `InferenceSurvivalWeibullRegr` | Survival | Likelihood | Not implemented | Weibull AFT |
| `InferenceSurvivalDepCensTransformRegr` | Survival | Likelihood | Not implemented | Dependent-censoring transform |
| `InferenceSurvivalKKRankRegrIVWC` | Survival | KK IVWC | Not implemented | Rank-regression IVWC |
| `InferenceSurvivalKKLWACoxIVWC` | Survival | KK IVWC likelihood | Not implemented | LWA Cox IVWC |
| `InferenceSurvivalKKLWACoxOneLik` | Survival | KK one-likelihood | Not implemented | LWA Cox one-likelihood |
| `InferenceSurvivalKKStratCoxIVWC` | Survival | KK IVWC likelihood | Not implemented | Stratified Cox IVWC |
| `InferenceSurvivalKKStratCoxOneLik` | Survival | KK one-likelihood | Not implemented | Stratified Cox one-likelihood |
| `InferenceSurvivalKKWeibullFrailtyIVWC` | Survival | KK IVWC likelihood | Not implemented | Weibull frailty IVWC |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | Survival | KK one-likelihood | Not implemented | Weibull frailty one-likelihood |
| `InferenceSurvivalKKClaytonCopulaIVWC` | Survival | KK IVWC likelihood | Not implemented | Clayton copula IVWC |
| `InferenceSurvivalKKClaytonCopulaOneLik` | Survival | KK one-likelihood | Not implemented | Clayton copula one-likelihood |

## Out of Scope for the Current Weighted Bayesian-Bootstrap Rollout

| Class | Family | Path type | Status | Notes |
|---|---|---|---|---|
| `InferenceIncidExactFisher` | Incidence | Exact | Out of scope | Exact-only path |
| `InferenceIncidenceExactBinomial` | Incidence | Exact | Out of scope | Exact-only path |
| `InferenceIncidExactZhang` | Incidence | Exact | Out of scope | Exact / inversion path |
| `InferenceCustomAsymp` | Custom | Extension | Out of scope | User extension path |
| `InferenceCustomRand` | Custom | Extension | Out of scope | User extension path |
| `InferenceCustomBoot` | Custom | Extension | Out of scope | User extension path |

## Summary

- Implemented concrete paths: `3`
- Not yet implemented concrete paths tracked here: `96`
- Out-of-scope concrete paths tracked here: `6`

Immediate next-wave candidates are likely:

1. Simple empirical estimators:
   `InferenceAllSimpleMeanDiff`, `InferenceAllSimpleMeanDiffPooledVar`,
   `InferenceIncidRiskDiff`, `InferenceSurvivalKMDiff`,
   `InferenceSurvivalRestrictedMeanDiff`
2. G-computation paths:
   `InferenceIncidGCompRiskDiff`, `InferenceIncidGCompRiskRatio`,
   `InferencePropGCompMeanDiff`, `InferenceOrdinalGCompMeanDiff`
3. GEE paths:
   `InferenceIncidKKGEE`, `InferenceCountPoissonKKGEE`,
   `InferenceCountPoissonUnivKKGEE`, `InferenceCountPoissonMultiKKGEE`,
   `InferenceOrdinalKKGEE`, `InferencePropKKGEE`
