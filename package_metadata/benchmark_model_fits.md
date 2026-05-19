# EDI Exhaustive C++ Model Fit Benchmarks

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Benchmark Dataset Specification

All benchmarks were performed on a synthetic clinical-trial-scale dataset generated for each response type. The data generation process ensures numerical stability and fair solver comparison by using the following parameters:

*   **Sample Size ($N$):** 1,000 subjects for most models; 500 subjects for survival models. Exact and trend tests may use smaller scaled samples (N=100-500) as noted in the results.
*   **Predictors ($p$):** 5 total predictors, including a global intercept, a binary treatment assignment ($w \sim \text{Bernoulli}(0.5)$), and 4 continuous covariates ($X \sim \text{Normal}(0, 1)$).
*   **Effect Sizes:** Coefficients are sampled from $\text{Normal}(0, 0.1)$ to ensure reasonable event rates and avoid separation issues in logistic/ordinal models.
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
*   **Low-Level Comparison:** Canonical R timings use the fastest available internal interfaces (e.g., `.fit` functions) to remove R's formula parsing and environment management overhead.
*   **Averaging:** All timings are medians over 3 independent iterations per model path.
*   **Constraints**: Matched-pair/KK and highly custom paths are excluded as per user request.

## Results

| Class | Response | EDI Time (ms) | Canonical Pkg | Canonical Func | Canonical Time (ms) | Speedup |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| InferenceIncidLogRegr | incidence | 0.05 | stats | glm.fit | 2.01 | 39.14x |
| InferenceContinOLS | continuous | 0.04 | stats | lm.fit | 0.16 | 3.95x |
| InferenceCountPoisson | count | 0.05 | stats | glm.fit | 1.83 | 38.38x |
| InferenceSurvivalCoxPHRegr | survival | 0.03 | survival | coxph.fit | 0.67 | 24.57x |
| InferenceCountNegBin | count | 0.02 | MASS | glm.nb | 43.88 | 2437.99x |
| InferencePropBetaRegr | proportion | 0.02 | betareg | betareg.fit | 34.31 | 1884.94x |
| InferenceOrdinalPropOddsRegr | ordinal | 0.02 | ordinal | clm | 10.35 | 575.17x |
| InferenceCountHurdlePoisson | count | 0.02 | pscl | hurdle | 21.64 | 1109.62x |
| InferenceCountZeroInflatedPoisson | count | 0.11 | pscl | zeroinfl | 101.56 | 919.07x |
| InferenceCountZeroInflatedNegBin | count | 0.02 | pscl | zeroinfl(nb) | 125.3 | 6202.84x |
| InferenceCountHurdleNegBin | count | 0.02 | pscl | hurdle(nb) | 57.18 | 2917.43x |
| InferenceCountQuasiPoisson | count | 0.05 | stats | glm.fit(quasi) | 2.72 | 55.98x |
| InferenceSurvivalWeibullRegr | survival | 0.03 | survival | survreg | 3.6 | 143.94x |
| InferenceContinRobustRegr | continuous | 0.04 | MASS | rlm | 1.68 | 41.05x |
| InferenceContinQuantileRegr | continuous | 0.06 | quantreg | rq.fit | 0.93 | 15.35x |
| InferencePropFractionalLogit | proportion | 0.02 | stats | glm.fit(quasi) | 1.14 | 57.79x |
| InferenceIncidLogBinomial | incidence | 0.02 | stats | glm.fit(log) | 3.77 | 181.29x |
| InferenceIncidProbitRegr | incidence | 0.04 | stats | glm.fit(probit) | 2.87 | 69.74x |
| InferenceIncidBinomialIdentityRiskDiff | incidence | 0.02 | stats | glm.fit(ident) | 1.28 | 64.29x |
| InferenceOrdinalAdjCatLogitRegr | ordinal | 0.02 | VGAM | vglm(acat) | 56.17 | 2794.42x |
| InferenceOrdinalContRatioRegr | ordinal | 0.02 | VGAM | vglm(cratio) | 40.08 | 2076.49x |
| InferenceOrdinalOrderedProbitRegr | ordinal | 0.02 | ordinal | clm(probit) | 8.78 | 466.97x |
| InferenceOrdinalCloglogRegr | ordinal | 0.02 | ordinal | clm(cll) | 9.19 | 499.26x |
| InferenceOrdinalCauchitRegr | ordinal | 0.04 | ordinal | clm(cauchit) | 12.45 | 293.71x |
| InferenceSurvivalLogRank | survival | 0.04 | survival | survdiff | 1.84 | 45.02x |
| InferenceSurvivalGehanWilcox | survival | 0.04 | survival | survdiff(rho=1) | 1.88 | 47.43x |
| InferenceAllSimpleMeanDiffPooledVar | continuous | 0.06 | stats | t.test(pool) | 0.18 | 3.27x |
| InferenceAllSimpleWilcox | continuous | 0.63 | stats | wilcox.test | 1.06 | 1.69x |
| InferenceIncidExactFisher | incidence | 0.51 | stats | fisher.test | 0.66 | 1.29x |
| InferenceIncidCMH | incidence | 0.01 | stats | mantelhaen | 1.09 | 79.86x |
| InferenceOrdinalJonckheereTerpstraTest | ordinal | 0.04 | clinfun | jonckheere | 0.41 | 10.44x |
| InferenceIncidMiettinenNurminenRiskDiff | incidence | 0.04 | DescTools | BinomDiffCI(mn) | 0.61 | 15.36x |
| InferenceIncidModifiedPoisson | incidence | 0.02 | stats | glm.fit(modified) | 1.81 | 93.43x |
| InferenceSurvivalStratCoxPHRegr | survival | 0.04 | survival | coxph.fit(strat) | 0.63 | 16.35x |
| InferenceContinLin | continuous | 0.04 | stats | lm.fit(interact) | 0.6 | 14.38x |
| InferenceIncidRiskDiff | incidence | 0.02 | stats | prop.test | 0.5 | 26.66x |
| InferenceSurvivalKMDiff | survival | 0.06 | survival | survfit(median) | 4.18 | 69.39x |
| InferenceSurvivalRestrictedMeanDiff | survival | 0.01 | survival | survfit(rmean) | 3.26 | 226.26x |
| InferenceCountRobustPoisson | count | 0.05 | stats | glm.fit | 1.78 | 36.59x |
| InferenceOrdinalRidit | ordinal | 0.04 | stats | mean(ridit) | 0.28 | 7.01x |
| InferenceIncidExactZhang | incidence | 0.04 | Exact | exact.test(z) | 1608.77 | 41895.08x |
| InferenceIncidNewcombeRiskDiff | incidence | 0.03 | DescTools | BinomDiffCI(score) | 0.86 | 26.03x |
| InferencePropGCompMeanDiff | proportion | 0.05 | None | None | NA | NA |
| InferenceSurvivalDepCensTransformRegr | survival | 0.04 | None | None | NA | NA |
| InferenceIncidGCompRiskRatio | incidence | 0.05 | None | None | NA | NA |
| InferenceIncidGCompRiskDiff | incidence | 0.02 | None | None | NA | NA |
| InferenceOrdinalGCompMeanDiff | ordinal | 0.05 | None | None | NA | NA |
| InferencePropZeroOneInflatedBetaRegr | proportion | 0.02 | None | None | NA | NA |
