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

- In-scope concrete bootstrap-capable classes with real weighted support: `31`
- In-scope concrete bootstrap-capable classes still unimplemented: `77`
- Exact/custom out-of-scope classes: `6`
- Non-IVWC unimplemented classes shown in the main tables below: `64`
- IVWC classes omitted from the main tables: `13`

## Crosscheck Against The Previous Checklist

The previous version of this file was stale.

Main discrepancies:

- It said only `3` concrete paths were implemented. The current source has `31`.
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
  `InferenceIncidMultiModifiedPoisson`,
  `InferenceOrdinalStereotypeLogitRegr`,
  `InferenceOrdinalUniStereotypeLogitRegr`,
  `InferenceOrdinalMultiStereotypeLogitRegr`,
  `InferenceOrdinalContRatioRegr`,
  `InferenceOrdinalUniContRatioRegr`,
  `InferenceOrdinalMultiContRatioRegr`,
  `InferenceOrdinalUniCumulProbitRegr`,
  `InferenceOrdinalMultiCumulProbitRegr`.
- It listed custom extension classes in the tracked table even though they should
  remain out of scope for the package’s built-in weighted rollout.

## Implemented

These classes currently have a real weighted implementation, either directly or
through a non-stub parent class.

Two additional implemented paths are omitted from the main table because they
are `IVWC`:

- `InferenceCountKKCPoissonIVWC`
- `InferenceCountKKHurdlePoissonIVWC`

### Implemented Non-IVWC Paths

| Class | Family | Provider | Notes |
|---|---|---|---|
| `InferenceAllSimpleMeanDiff` | All / continuous | `InferenceAllSimpleMeanDiff` | Weighted empirical mean difference |
| `InferenceAllSimpleMeanDiffPooledVar` | All / continuous | `InferenceAllSimpleMeanDiff` | Inherits weighted mean-difference estimate |
| `InferenceContinLin` | Continuous | `InferenceContinLin` | Weighted `lm.wfit` on Lin design |
| `InferenceContinOLS` | Continuous | `InferenceContinOLS` | Weighted `lm.wfit` |
| `InferenceCountKKCPoissonOneLik` | Count KK one-likelihood | `InferenceAbstractKKPoissonCPoissonOneLik` | Weighted one-likelihood path |
| `InferenceCountKKHurdlePoissonOneLik` | Count KK one-likelihood | `InferenceAbstractKKHurdlePoissonOneLik` | Weighted one-likelihood path |
| `InferenceCountPoisson` | Count | `InferenceCountPoisson` | First-wave weighted Poisson |
| `InferenceCountPoissonKKGEE` | Count KK GEE | `InferenceCountPoissonKKGEE` | Native weighted KK GEE |
| `InferenceCountPoissonMultiKKGEE` | Count KK GEE | `InferenceCountPoissonKKGEE` | Inherits weighted KK GEE implementation |
| `InferenceCountPoissonUnivKKGEE` | Count KK GEE | `InferenceCountPoissonKKGEE` | Inherits weighted KK GEE implementation |
| `InferenceIncidCMH` | Incidence | `InferenceAllSimpleMeanDiff` | Weighted empirical contrast path |
| `InferenceIncidExtendedRobins` | Incidence | `InferenceAllSimpleMeanDiff` | Weighted empirical contrast path |
| `InferenceIncidGCompRiskDiff` | Incidence g-comp | `InferenceIncidGCompAbstract` | Weighted g-computation |
| `InferenceIncidGCompRiskRatio` | Incidence g-comp | `InferenceIncidGCompAbstract` | Weighted g-computation |
| `InferenceIncidKKGEE` | Incidence KK GEE | `InferenceIncidKKGEE` | Native weighted KK GEE |
| `InferenceIncidLogRegr` | Incidence | `InferenceIncidLogRegr` | First-wave weighted logistic |
| `InferenceIncidModifiedPoisson` | Incidence | `InferenceIncidModifiedPoisson` | First-wave weighted modified Poisson |
| `InferenceIncidMultiModifiedPoisson` | Incidence | `InferenceIncidModifiedPoisson` | Alias inherits weighted modified Poisson |
| `InferenceIncidRiskDiff` | Incidence | `InferenceIncidRiskDiff` | Weighted identity-scale risk difference |
| `InferenceIncidenceWald` | Incidence | `InferenceAllSimpleMeanDiff` | Weighted empirical contrast path |
| `InferenceOrdinalGCompMeanDiff` | Ordinal g-comp | `InferenceOrdinalGCompMeanDiff` | Weighted ordinal g-computation |
| `InferenceOrdinalKKGEE` | Ordinal KK GEE | `InferenceOrdinalKKGEE` | Weighted ordinal surrogate path |
| `InferenceOrdinalPropOddsRegr` | Ordinal | `InferenceOrdinalPropOddsRegr` | Weighted proportional odds |
| `InferencePropGCompMeanDiff` | Proportion g-comp | `InferencePropGCompAbstract` | Weighted g-computation |
| `InferencePropKKGEE` | Proportion KK GEE | `InferencePropKKGEE` | Native weighted KK GEE |
| `InferenceSurvivalGehanWilcox` | Survival | `InferenceSurvivalGehanWilcox` | Weighted Peto-Prentice score estimate |
| `InferenceSurvivalKMDiff` | Survival | `InferenceSurvivalKMDiff` | Weighted KM median difference |
| `InferenceSurvivalLogRank` | Survival | `InferenceSurvivalLogRank` | Weighted log-rank score estimate |
| `InferenceSurvivalRestrictedMeanDiff` | Survival | `InferenceSurvivalRestrictedMeanDiff` | Weighted RMST difference |

## Still Unimplemented

These are concrete bootstrap-capable classes that still fall through to the
Bayesian-bootstrap stub and therefore do not yet have a real weighted
re-estimation path.

### Near-Term / Moderate Difficulty

These are the most plausible next implementations if we continue expanding
coverage without adding major new infrastructure.

| Class | Family | Why It Is Not Done Yet |
|---|---|---|
| `InferenceContinRobustRegr` | Continuous | Needs a settled policy for observation-weighted robust regression |
| `InferenceContinQuantileRegr` | Continuous | Needs weighted quantile-regression implementation |
| `InferenceCountNegBin` | Count | Needs weighted negative-binomial refit |
| `InferenceCountQuasiPoisson` | Count | Needs weighted quasi-likelihood path |
| `InferenceCountRobustPoisson` | Count | Needs weighted robust/sandwich recomputation |
| `InferenceCountHurdlePoisson` | Count | Two-part weighted mixture refit needed |
| `InferenceCountHurdleNegBin` | Count | Two-part weighted mixture refit needed |
| `InferenceCountZeroInflatedPoisson` | Count | Zero-inflated weighted refit needed |
| `InferenceCountZeroInflatedNegBin` | Count | Zero-inflated weighted refit needed |
| `InferenceIncidBinomialIdentityRiskDiff` | Incidence | Weighted constrained binomial path not yet wired |
| `InferenceIncidLogBinomial` | Incidence | Weighted constrained log-binomial path not yet wired |
| `InferenceIncidProbitRegr` | Incidence | No weighted probit backend in current native stack |
| `InferenceIncidNewcombeRiskDiff` | Incidence | Standalone non-likelihood formula path not yet adapted |
| `InferenceIncidMiettinenNurminenRiskDiff` | Incidence | Standalone non-likelihood formula path not yet adapted |
| `InferencePropFractionalLogit` | Proportion | Needs weighted fractional-logit refit |
| `InferencePropBetaRegr` | Proportion | Needs weighted beta-regression refit |
| `InferencePropZeroOneInflatedBetaRegr` | Proportion | Mixture-model weighted refit needed |

### Ordinal Likelihood Gaps

These remain unimplemented largely because the current native weighted ordinal
backend only covers the proportional-odds logit-link path.

| Class | Family | Why It Is Not Done Yet |
|---|---|---|
| `InferenceOrdinalAdjCatLogitRegr` | Ordinal | No native weighted adjacent-category backend |
| `InferenceOrdinalCloglogRegr` | Ordinal | No native weighted cloglog ordinal backend |
| `InferenceOrdinalCauchitRegr` | Ordinal | No native weighted cauchit ordinal backend |
| `InferenceOrdinalOrderedProbitRegr` | Ordinal | No native weighted ordered-probit backend |
| `InferenceOrdinalStereotypeLogitRegr` | Ordinal | No weighted stereotype-logit implementation |
| `InferenceOrdinalUniStereotypeLogitRegr` | Ordinal | Inherits unimplemented stereotype-logit path |
| `InferenceOrdinalMultiStereotypeLogitRegr` | Ordinal | Inherits unimplemented stereotype-logit path |
| `InferenceOrdinalContRatioRegr` | Ordinal | No weighted continuation-ratio backend |
| `InferenceOrdinalUniContRatioRegr` | Ordinal | Inherits unimplemented continuation-ratio path |
| `InferenceOrdinalMultiContRatioRegr` | Ordinal | Inherits unimplemented continuation-ratio path |
| `InferenceOrdinalUniCumulProbitRegr` | Ordinal | Inherits unimplemented cumulative-probit path |
| `InferenceOrdinalMultiCumulProbitRegr` | Ordinal | Inherits unimplemented cumulative-probit path |
| `InferenceOrdinalRidit` | Ordinal | Non-likelihood weighted path not yet defined |
| `InferenceOrdinalPairedSignTest` | Ordinal | Paired sign-test weighted analogue not yet defined |

### KK One-Likelihood / GLMM / Compound Paths

These are bootstrap-capable, but weighted observation support likely requires
new numerical work or a more careful definition of the weighted empirical
estimand.

| Class | Family |
|---|---|
| `InferenceContinKKGLMM` | Continuous KK GLMM |
| `InferenceContinKKOLSOneLik` | Continuous KK one-likelihood |
| `InferenceContinKKQuantileRegrOneLik` | Continuous KK one-likelihood |
| `InferenceContinKKRobustRegrOneLik` | Continuous KK one-likelihood |
| `InferenceCountCompositeLikelihood` | Count composite likelihood |
| `InferenceCountKKGLMM` | Count KK GLMM |
| `InferenceIncidKKClogitOneLik` | Incidence KK one-likelihood |
| `InferenceIncidKKClogitPlusGLMMOneLik` | Incidence KK hybrid one-likelihood |
| `InferenceIncidKKGCompRiskDiff` | Incidence KK g-comp |
| `InferenceIncidKKGCompRiskRatio` | Incidence KK g-comp |
| `InferenceIncidKKGLMM` | Incidence KK GLMM |
| `InferenceIncidKKModifiedPoisson` | Incidence KK marginal |
| `InferenceIncidKKNewcombeRiskDiff` | Incidence KK non-likelihood |
| `InferenceOrdinalKKCondAdjCatLogitRegr` | Ordinal KK likelihood |
| `InferenceOrdinalKKGLMM` | Ordinal KK GLMM |
| `InferencePropKKGLMM` | Proportion KK GLMM |
| `InferencePropKKQuantileRegrOneLik` | Proportion KK one-likelihood |

### Survival Regression / Survival KK Deferred Paths

These are the heaviest remaining paths and are the clearest candidates to defer
until after simpler GLM and empirical-function coverage is complete.

| Class | Family |
|---|---|
| `InferenceSurvivalCoxPHRegr` | Survival Cox PH |
| `InferenceSurvivalStratCoxPHRegr` | Survival stratified Cox PH |
| `InferenceSurvivalWeibullRegr` | Survival Weibull |
| `InferenceSurvivalDepCensTransformRegr` | Survival dependent-censoring transform |
| `InferenceSurvivalKKClaytonCopulaOneLik` | Survival KK copula one-likelihood |
| `InferenceSurvivalKKLWACoxOneLik` | Survival KK LWA Cox one-likelihood |
| `InferenceSurvivalKKStratCoxOneLik` | Survival KK stratified Cox one-likelihood |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | Survival KK Weibull frailty one-likelihood |

## IVWC Omitted From The Main Tables

These bootstrap-capable classes are still unimplemented, but are omitted from
the main tables above per prior preference to leave IVWC paths out of the main
inventory.

- `InferenceAllKKMeanDiffIVWC`
- `InferenceAllKKWilcoxIVWC`
- `InferenceContinKKOLSIVWC`
- `InferenceContinKKQuantileRegrIVWC`
- `InferenceContinKKRobustRegrIVWC`
- `InferenceIncidKKClogitIVWC`
- `InferenceIncidKKClogitPlusGLMMIVWC`
- `InferencePropKKQuantileRegrIVWC`
- `InferenceSurvivalKKClaytonCopulaIVWC`
- `InferenceSurvivalKKLWACoxIVWC`
- `InferenceSurvivalKKRankRegrIVWC`
- `InferenceSurvivalKKStratCoxIVWC`
- `InferenceSurvivalKKWeibullFrailtyIVWC`

## Out Of Scope

These are not part of the package’s built-in weighted Bayesian-bootstrap rollout.

- `InferenceIncidExactFisher`
- `InferenceIncidenceExactBinomial`
- `InferenceIncidExactZhang`
- `InferenceCustomAsymp`
- `InferenceCustomRand`
- `InferenceCustomBoot`
