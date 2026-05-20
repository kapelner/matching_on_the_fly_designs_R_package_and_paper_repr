# EDI Exhaustive C++ Model Fit Benchmarks

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Benchmark Dataset Specification

All benchmarks were performed on a synthetic clinical-trial-scale dataset generated for each response type. The data generation process ensures numerical stability and fair solver comparison by using the following parameters:

*   **Sample Size ($N$):** 1,000 subjects for most models; 500 subjects for survival models. Exact and trend tests may use smaller scaled samples (N=100-500) as noted in the results.
*   **Predictors ($p$):** 5 total predictors, including a global intercept, a balanced binary treatment assignment from fixed `iBCRD`, and 4 continuous covariates ($X \sim \text{Normal}(0, 1)$).
*   **Effect Sizes:** Coefficients are sampled from $\text{Normal}(0, 0.1)$ to ensure reasonable event rates and avoid separation issues in logistic/ordinal models.
*   **EDI Design Template:** EDI benchmark objects are instantiated on a fixed `iBCRD` design.
*   **Response Generation:**
    *   **Continuous:** Linear model with additive $\text{Normal}(0, 0.5)$ noise.
    *   **Incidence:** Binary outcomes via a Logistic link.
    *   **Count:** Integer outcomes via Poisson or Negative Binomial distributions with an exponential link.
    *   **Proportion:** Continuous outcomes in $(0, 1)$ via a Beta distribution with a logit link.
    *   **Survival:** Exponentially distributed event times with approximately 20% random censoring.
    *   **Ordinal:** 4-level categorical outcomes generated from a Proportional Odds model.

## Methodology

*   **Pure Solver Timing:** Results reflect the time taken for the core numerical optimization. We exclude R6 object instantiation, design matrix construction, and standard error estimation (which often uses different R-side matrix inversion logic) to isolate the efficiency of the underlying C++ backends.
*   **Smart Cold Starts:** EDI models were initialized with `smart_cold_start = TRUE`, utilizing package-optimized heuristic starting values.
*   **Randomization Design:** EDI timings in this table correspond to `iBCRD` design objects.
*   **Low-Level Comparison:** Canonical R timings use the fastest available internal interfaces (e.g., `.fit` functions) to remove R's formula parsing and environment management overhead.
*   **Averaging:** All timings are medians over 3 warmed runs measured with adaptive batched `system.time`; paths below 0.01 ms use `microbenchmark(times = 500)` instead.
*   **Constraints**: Matched-pair/KK and highly custom paths are excluded as per user request.

## Results

| Class | Response | EDI Time (ms) | Canonical Pkg | Canonical Func | Canonical Time (ms) | Speedup |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| InferenceAllSimpleMeanDiffPooledVar | continuous | 0.05 | stats | t.test(pool) | 0.18 | 3.96x |
| InferenceAllSimpleWilcox | continuous | 0.51 | stats | wilcox.test | 1.11 | 2.17x |
| InferenceContinLin | continuous | 1.00 | stats | lm.fit(interact) | 0.73 | 0.73x |
| InferenceContinOLS | continuous | 0.50 | stats | lm.fit | 0.13 | 0.25x |
| InferenceContinQuantileRegr | continuous | 2.00 | quantreg | rq.fit | 0.93 | 0.47x |
| InferenceContinRobustRegr | continuous | 0.50 | MASS | rlm | 2.00 | 4x |
| InferenceIncidBinomialIdentityRiskDiff | incidence | 1.50 | stats | glm.fit(ident) | 1.73 | 1.15x |
| InferenceIncidExactFisher | incidence | 0.67 | stats | fisher.test | 0.88 | 1.32x |
| InferenceIncidExactZhang | incidence | 0.04 | Exact | exact.test(z) | 1980.00 | 47826.09x |
| InferenceIncidGCompRiskDiff | incidence | 1.05 | stats | glm.fit+gcomp(RD) | 1.87 | 1.78x |
| InferenceIncidGCompRiskRatio | incidence | 1.00 | stats | glm.fit+gcomp(RR) | 1.37 | 1.37x |
| InferenceIncidLogBinomial | incidence | 1.00 | stats | glm.fit(log) | 1.68 | 1.68x |
| InferenceIncidLogRegr | incidence | 1.00 | stats | glm.fit | 2.59 | 2.59x |
| InferenceIncidMiettinenNurminenRiskDiff | incidence | 1.00 | DescTools | BinomDiffCI(mn) | 1.33 | 1.33x |
| InferenceIncidModifiedPoisson | incidence | 1.50 | stats | glm.fit(modified) | 2.80 | 1.87x |
| InferenceIncidNewcombeRiskDiff | incidence | 0.04 | DescTools | BinomDiffCI(score) | 1.22 | 33.26x |
| InferenceIncidProbitRegr | incidence | 2.50 | stats | glm.fit(probit) | 1.96 | 0.78x |
| InferenceIncidRiskDiff | incidence | 0.50 | stats | prop.test | 0.90 | 1.81x |
| InferenceCountHurdleNegBin | count | 4.00 | pscl | hurdle(nb) | 44.00 | 11x |
| InferenceCountHurdlePoisson | count | 3.00 | pscl | hurdle | 23.00 | 7.67x |
| InferenceCountNegBin | count | 1.50 | MASS | glm.nb | 41.00 | 27.33x |
| InferenceCountPoisson | count | 1.50 | stats | glm.fit | 2.00 | 1.33x |
| InferenceCountQuasiPoisson | count | 0.50 | stats | glm.fit(quasi) | 2.19 | 4.37x |
| InferenceCountRobustPoisson | count | 2.00 | stats | glm.fit | 2.26 | 1.13x |
| InferenceCountZeroInflatedNegBin | count | 2.84 | pscl | zeroinfl(nb) | 171.00 | 60.17x |
| InferenceCountZeroInflatedPoisson | count | 2.26 | pscl | zeroinfl | 92.00 | 40.65x |
| InferencePropBetaRegr | proportion | 8.00 | betareg | betareg.fit | 37.50 | 4.69x |
| InferencePropFractionalLogit | proportion | 2.50 | stats | glm.fit(quasi) | 1.23 | 0.49x |
| InferencePropGCompMeanDiff | proportion | 1.00 | stats | glm.fit(quasi)+gcomp | 1.41 | 1.41x |
| InferencePropZeroOneInflatedBetaRegr | proportion | 3.50 | None | None | NA | NA |
| InferenceSurvivalCoxPHRegr | survival | 44.00 | survival | coxph.fit | 0.59 | 0.01x |
| InferenceSurvivalDepCensTransformRegr | survival | 9.00 | None | None | NA | NA |
| InferenceSurvivalGehanWilcox | survival | 5.50 | survival | survdiff(rho=1) | 2.15 | 0.39x |
| InferenceSurvivalKMDiff | survival | 0.05 | survival | survfit(median) | 4.21 | 91.22x |
| InferenceSurvivalLogRank | survival | 0.25 | survival | survdiff | 2.50 | 10x |
| InferenceSurvivalRestrictedMeanDiff | survival | 0.05 | survival | survfit(rmean) | 2.87 | 55.37x |
| InferenceSurvivalStratCoxPHRegr | survival | 82.00 | survival | coxph.fit(strat) | 1.05 | 0.01x |
| InferenceSurvivalWeibullRegr | survival | 1.00 | survival | survreg | 3.87 | 3.87x |
| InferenceOrdinalAdjCatLogitRegr | ordinal | 3.67 | VGAM | vglm(acat) | 50.00 | 13.64x |
| InferenceOrdinalCauchitRegr | ordinal | 7.00 | ordinal | clm(cauchit) | 16.25 | 2.32x |
| InferenceOrdinalCloglogRegr | ordinal | 6.00 | ordinal | clm(cll) | 11.60 | 1.93x |
| InferenceOrdinalContRatioRegr | ordinal | 3.00 | VGAM | vglm(cratio) | 44.00 | 14.67x |
| InferenceOrdinalGCompMeanDiff | ordinal | 6.00 | ordinal | clm+gcomp | 23.67 | 3.94x |
| InferenceOrdinalJonckheereTerpstraTest | ordinal | 39.00 | clinfun | jonckheere | 1.32 | 0.03x |
| InferenceOrdinalOrderedProbitRegr | ordinal | 4.50 | ordinal | clm(probit) | 9.33 | 2.07x |
| InferenceOrdinalPropOddsRegr | ordinal | 5.50 | ordinal | clm | 11.25 | 2.05x |
| InferenceOrdinalRidit | ordinal | 0.06 | stats | mean(ridit) | 0.32 | 5.04x |
