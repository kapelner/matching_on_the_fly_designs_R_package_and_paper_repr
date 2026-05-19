# EDI Exhaustive Model Fit Benchmarks

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.
Timings represent pure solver execution (excluding R6 object instantiation overhead) with `smart_cold_start = TRUE` enabled for EDI.

| Class | Response | EDI Time (ms) | Canonical Pkg | Canonical Func | Canonical Time (ms) | Speedup |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| InferenceIncidLogRegr | incidence | 29.88 | stats | glm.fit | 2.5 | 0.08x |
| InferenceContinOLS | continuous | 0.85 | stats | lm.fit | 0.32 | 0.38x |
| InferenceCountPoisson | count | 5.77 | stats | glm.fit | 2.22 | 0.38x |
| InferenceSurvivalCoxPHRegr | survival | 0.04 | survival | coxph.fit | 1.04 | 23.84x |
| InferenceCountNegBin | count | 0.04 | MASS | glm.nb | 64.21 | 1472.7x |
| InferencePropBetaRegr | proportion | 0.04 | betareg | betareg.fit | 87.93 | 2134.21x |
| InferenceOrdinalPropOddsRegr | ordinal | 0.07 | ordinal | clm | 19.29 | 295.92x |
| InferenceCountHurdlePoisson | count | 0.1 | pscl | hurdle | 34.43 | 339.84x |
| InferenceCountZeroInflatedPoisson | count | 0.05 | pscl | zeroinfl | 143.6 | 3029.47x |
| InferenceCountZeroInflatedNegBin | count | 0.07 | pscl | zeroinfl(nb) | 197.59 | 2949.05x |
| InferenceCountHurdleNegBin | count | 0.08 | pscl | hurdle(nb) | 79.75 | 1026.38x |
| InferenceCountQuasiPoisson | count | 9.13 | stats | glm.fit(quasi) | 2.26 | 0.25x |
| InferenceSurvivalWeibullRegr | survival | 0.04 | survival | survreg | 6.32 | 153.79x |
| InferenceContinRobustRegr | continuous | 16.22 | MASS | rlm | 2.76 | 0.17x |
| InferenceContinQuantileRegr | continuous | 7.09 | quantreg | rq.fit | 1.3 | 0.18x |
| InferencePropFractionalLogit | proportion | 0.04 | stats | glm.fit(quasi) | 2.01 | 45.73x |
| InferenceIncidLogBinomial | incidence | 0.05 | stats | glm.fit(log) | 3.54 | 71.03x |
| InferenceIncidProbitRegr | incidence | 0.05 | stats | glm.fit(probit) | 2.45 | 52.08x |
| InferenceIncidBinomialIdentityRiskDiff | incidence | 0.06 | stats | glm.fit(ident) | 2.88 | 49.77x |
| InferenceOrdinalAdjCatLogitRegr | ordinal | 0.07 | VGAM | vglm(acat) | 106.79 | 1593.93x |
| InferenceOrdinalContRatioRegr | ordinal | 0.07 | VGAM | vglm(cratio) | 78.07 | 1144.67x |
| InferenceOrdinalOrderedProbitRegr | ordinal | 0.38 | ordinal | clm(probit) | 15.83 | 41.78x |
| InferenceOrdinalCloglogRegr | ordinal | 0.05 | ordinal | clm(cll) | 17.42 | 336.97x |
| InferenceOrdinalCauchitRegr | ordinal | 0.05 | ordinal | clm(cauchit) | 20.63 | 382.79x |
| InferenceSurvivalLogRank | survival | 4.16 | survival | survdiff | 3.04 | 0.73x |
| InferenceSurvivalGehanWilcox | survival | 14.48 | survival | survdiff(rho=1) | 3.59 | 0.25x |
| InferenceAllSimpleMeanDiff | continuous | 0.43 | base | mean diff | 0.06 | 0.15x |
| InferenceAllSimpleMeanDiffPooledVar | continuous | 0.04 | stats | t.test(pool) | 1.32 | 30.1x |
| InferenceAllSimpleWilcox | continuous | 0.92 | stats | wilcox.test | 1.6 | 1.74x |
| InferenceIncidExactFisher | incidence | 0.74 | stats | fisher.test | 0.9 | 1.22x |
| InferenceIncidCMH | incidence | 0.02 | stats | mantelhaen | 1.28 | 63.23x |
| InferenceOrdinalPairedSignTest | ordinal | NA | stats | binom.test | 0.65 | NA |
| InferenceIncidModifiedPoisson | incidence | 0.04 | None | None | NA | NA |
| InferenceIncidenceWald | incidence | 0.02 | None | None | NA | NA |
| InferencePropGCompMeanDiff | proportion | 8.91 | None | None | NA | NA |
| InferenceIncidExtendedRobins | incidence | NA | None | None | NA | NA |
| InferenceIncidMiettinenNurminenRiskDiff | incidence | 0.27 | None | None | NA | NA |
| InferenceSurvivalStratCoxPHRegr | survival | 21.97 | None | None | NA | NA |
| InferenceContinLin | continuous | 14.94 | None | None | NA | NA |
| InferenceOrdinalJonckheereTerpstraTest | ordinal | 28.96 | None | None | NA | NA |
| InferenceIncidRiskDiff | incidence | 0.05 | None | None | NA | NA |
| InferenceSurvivalDepCensTransformRegr | survival | 0.2 | None | None | NA | NA |
| InferenceSurvivalKMDiff | survival | 0.1 | None | None | NA | NA |
| InferenceCountRobustPoisson | count | 0.13 | None | None | NA | NA |
| InferenceIncidNewcombeRiskDiff | incidence | 0.32 | None | None | NA | NA |
| InferenceIncidGCompRiskRatio | incidence | 13.16 | None | None | NA | NA |
| InferenceIncidExactBinomial | incidence | NA | None | None | NA | NA |
| InferenceIncidGCompRiskDiff | incidence | 0.09 | None | None | NA | NA |
| InferenceIncidExactZhang | incidence | 0.42 | None | None | NA | NA |
| InferenceOrdinalGCompMeanDiff | ordinal | 28.66 | None | None | NA | NA |
| InferencePropZeroOneInflatedBetaRegr | proportion | 0.06 | None | None | NA | NA |
| InferenceSurvivalRestrictedMeanDiff | survival | 0.03 | None | None | NA | NA |
| InferenceOrdinalRidit | ordinal | 0.31 | None | None | NA | NA |
