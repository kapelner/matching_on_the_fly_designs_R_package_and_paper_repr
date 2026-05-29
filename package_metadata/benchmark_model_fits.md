# EDI Exhaustive C++ Model Fit Benchmarks

_Generated: 2026-05-28 22:51:06 EDT_

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Compilation Context

These rows are read from build metadata compiled into the loaded `EDI` shared object via `edi_build_info_cpp()`.

**Compilation warning:** EDI model-fit timings are sensitive to the compiler flags used to build the loaded `EDI.so`. If EDI is compiled without the proper optimized flags, or with flags that are known to degrade these kernels such as problematic LTO builds, the benchmark can show substantial performance regressions that reflect the binary build rather than the modeling algorithms.

*   **EDI shared object:** `/home/kapelner/R/x86_64-pc-linux-gnu-library/4.7/EDI/libs/EDI.so`
*   **EDI shared object mtime:** `2026-05-28 22:22:44`
*   **Capture method:** `configure-generated header compiled into EDI.so`
*   **Build timestamp:** `2026-05-28 22:13:00 EDT`
*   **Build host:** `LAPTOP-J2T9TGGB`
*   **R version at build:** `R Under development (unstable) (2026-04-23 r89955) -- "Unsuffered Consequences"`
*   **R `CXX20` at build:** `g++ -O3 -march=native -flto -fno-math-errno`
*   **R `CXX20STD` at build:** `-std=gnu++20`
*   **R `CXX20FLAGS` at build:** `-g -O2`
*   **R `SHLIB_OPENMP_CXXFLAGS` at build:** `unavailable`
*   **Build env at build:** `EDI_PORTABLE=0`, `EDI_DISABLE_VECTORIZATION=0`, `EDI_NATIVE_SPEED=1`, `EDI_NATIVE_LTO=0`
*   **Package `PKG_CPPFLAGS` at build:** `-I../inst/include`
*   **Package `PKG_CXXFLAGS` at build:** `$(SHLIB_OPENMP_CXXFLAGS) -DNDEBUG -DEIGEN_NO_DEBUG -Wno-ignored-attributes -march=native -mtune=native -fno-lto; override CXXFLAGS+=-O3`
*   **Package `PKG_LIBS` at build:** `$(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -ltbb12 -fstack-protector`
*   **Compiler reported by binary:** `15.2.0`
*   **Compiler optimization macro enabled:** `TRUE`
*   **Compiler fast-math macro enabled:** `FALSE`
*   **Eigen vectorization disabled macro enabled:** `FALSE`

## Benchmark Dataset Specification

All benchmarks were performed on a synthetic clinical-trial-scale dataset generated for each response type. The data generation process ensures numerical stability and fair solver comparison by using the following parameters:

*   **Sample Size ($N$):** 1,000 subjects for most models; 500 subjects for survival models. Exact and trend tests may use smaller scaled samples (N=100-500) as noted in the results.
*   **Predictors ($p$):** 5 total predictors, including a global intercept, a balanced binary treatment assignment from fixed `iBCRD`, and 4 continuous covariates ($X \sim \text{Normal}(0, 1)$).
*   **Effect Sizes:** Covariate coefficients are sampled from $\text{Normal}(0, 0.5)$. The treatment coefficient is set to 0.5 in the linear predictor so the benchmarked treatment effect is meaningfully separated from zero.
*   **EDI Design Template:** EDI benchmark objects are instantiated on a fixed `iBCRD` design.
*   **Response Generation:**
    *   **Continuous:** Linear model with additive $\text{Normal}(0, 0.5)$ noise.
    *   **Incidence:** Binary outcomes via a Logistic link.
    *   **Count:** Integer outcomes via Poisson or Negative Binomial distributions with an exponential link.
    *   **Proportion:** Continuous outcomes in $(0, 1)$ via a Beta distribution with a logit link.
    *   **Survival:** Exponentially distributed event times with approximately 20% random censoring.
    *   **Ordinal:** 3-level categorical outcomes generated from the same ordinal construction used elsewhere in the benchmark suite.
*   **Stratified Cox Exception:** For `InferenceSurvivalStratCoxPHRegr`, the benchmark injects low-cardinality covariates before outcome generation so the row exercises a genuinely stratified Cox fit rather than the unstratified fallback.

## Methodology

*   **Bare Metal EDI Timing:** EDI rows call the exported C++ functions directly (e.g., `fast_logistic_regression_cpp`, `fast_ordinal_regression_cpp`) with all design matrices and fixed inputs pre-built outside the timed region. There is no R6 object instantiation, no cached state management, no warm start storage, and no standard error computation in the timed region — only the raw numerical solver.
*   **Apples-to-Apples Canonical Timing:** Canonical R timings likewise call the lowest-level publicly exposed interfaces (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) with pre-built design matrices. If a canonical package exposes no low-level function, the formula-based API is used instead.
*   **Low-Level Comparison:** Both EDI and canonical timings are measured on pre-built numeric matrices, removing formula parsing, model-frame construction, and R6/S3/S4 dispatch overhead from the timed region wherever the API permits.
*   **Limitation:** Some canonical comparators only expose formula-based APIs. Those rows remain included but their canonical timings carry formula/model-frame overhead not present in the EDI bare-metal timing.
*   **Averaging:** All timings are medians over 10 cold estimate-only timing samples measured with adaptive batched `system.time`; paths below 0.01 ms use `microbenchmark(times = 2000)` instead.
*   **Timing P-Value:** `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.
*   **Row Highlighting:** Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons.
*   **Constraints**: Matched-pair/KK and highly custom paths are excluded as per user request.

## Results

<table>
  <thead>
    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.82</td><td>stats</td><td>HL median pairwise diff</td><td>4.05</td><td>4.93x</td><td>8.59e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.02</td><td>stats</td><td>lm.fit</td><td>0.20</td><td>9.45x</td><td>3.4e-10</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>3.00</td><td>quantreg</td><td>rq.fit</td><td>2.50</td><td>0.83x</td><td>0.0602</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.26</td><td>MASS</td><td>rlm(MM)</td><td>56.67</td><td>216.05x</td><td>5.33e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.13</td><td>stats</td><td>glm.fit(ident)</td><td>24.25</td><td>187.44x</td><td>1.85e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.22</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>2.97</td><td>13.74x</td><td>1.5e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.25</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>3.16</td><td>12.82x</td><td>4.76e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>10.94</td><td>stats</td><td>glm.fit(log)</td><td>19.04</td><td>1.74x</td><td>1.82e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.17</td><td>stats</td><td>glm.fit</td><td>2.35</td><td>13.61x</td><td>1.41e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.27</td><td>stats</td><td>glm.fit(modified)</td><td>3.58</td><td>13.3x</td><td>1.84e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.79</td><td>stats</td><td>glm.fit(probit)</td><td>4.14</td><td>5.25x</td><td>7.04e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.02</td><td>stats</td><td>lm.fit(LPM)</td><td>0.26</td><td>12.89x</td><td>8.92e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>13.27</td><td>pscl</td><td>hurdle(nb)</td><td>84.50</td><td>6.37x</td><td>5.99e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>2.69</td><td>pscl</td><td>hurdle</td><td>35.80</td><td>13.33x</td><td>7.74e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>8.40</td><td>MASS</td><td>glm.nb</td><td>77.67</td><td>9.25x</td><td>6.54e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.21</td><td>stats</td><td>glm.fit</td><td>3.34</td><td>16.08x</td><td>1.99e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.22</td><td>stats</td><td>glm.fit(quasi)</td><td>3.66</td><td>17.04x</td><td>3.22e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.22</td><td>stats</td><td>glm.fit</td><td>3.26</td><td>14.82x</td><td>1.13e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>3.12</td><td>pscl</td><td>zeroinfl(nb)</td><td>280.00</td><td>89.72x</td><td>3.92e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>7.31</td><td>pscl</td><td>zeroinfl</td><td>118.75</td><td>16.24x</td><td>5.04e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>7.89</td><td>betareg</td><td>betareg.fit</td><td>52.13</td><td>6.61x</td><td>2.55e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.18</td><td>stats</td><td>glm.fit(quasi)</td><td>2.75</td><td>15.42x</td><td>1.1e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.23</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>2.92</td><td>12.76x</td><td>1.28e-09</td><td>***</td></tr>
    <tr><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.50</td><td>survival</td><td>coxph.fit(breslow)</td><td>1.10</td><td>2.2x</td><td>0.1</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.03</td><td>survival</td><td>survfit(median)</td><td>7.10</td><td>253.97x</td><td>9.2e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.03</td><td>survival</td><td>survdiff</td><td>3.61</td><td>110.36x</td><td>1.64e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.02</td><td>survival</td><td>survfit(rmean)</td><td>4.84</td><td>195x</td><td>5.99e-06</td><td>***</td></tr>
    <tr><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>1.00</td><td>survival</td><td>coxph.fit(strat)</td><td>0.68</td><td>0.68x</td><td>5.73e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.09</td><td>survival</td><td>survreg</td><td>6.63</td><td>74.96x</td><td>4.05e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.83</td><td>VGAM</td><td>vglm(acat)</td><td>49.88</td><td>60.13x</td><td>2.48e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>1.31</td><td>ordinal</td><td>clm(cauchit)</td><td>17.58</td><td>13.41x</td><td>7.58e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>1.25</td><td>ordinal</td><td>clm(cloglog)</td><td>16.80</td><td>13.4x</td><td>5.03e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.54</td><td>VGAM</td><td>vglm(cratio)</td><td>34.60</td><td>64.09x</td><td>4.64e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>1.32</td><td>ordinal</td><td>clm+gcomp</td><td>30.56</td><td>23.15x</td><td>5.84e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>1.26</td><td>ordinal</td><td>clm(probit)</td><td>13.88</td><td>10.98x</td><td>1.66e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>1.10</td><td>ordinal</td><td>clm</td><td>12.71</td><td>11.58x</td><td>1.3e-08</td><td>***</td></tr>
  </tbody>
</table>



## Wald Test Performance (Full Inference)

This table compares the performance of **Full Inference** (Model Fit + Standard Error calculation + P-value derivation).
Unlike the point-estimation table above, these results include the computational cost of the variance-covariance matrix (Hessian or Fisher Information) and the Wald test statistic calculation.
All paths (EDI and Canonical) use a reduced sample size ($N=200$) for this full-inference benchmark to ensure iterative stability.
EDI timings in this table correspond to fixed `iBCRD` design objects.
**Stratified Cox Exception**: For `InferenceSurvivalStratCoxPHRegr`, the benchmark injects low-cardinality covariates before outcome generation so the row exercises a genuinely stratified Cox fit rather than the unstratified fallback.
EDI regression models (Logistic, Poisson) are benchmarked using the **IRLS** optimizer for these Wald tests.
**Note on Coverage**: `InferenceIncidCMH` is retained for coverage, but under `iBCRD` its EDI asymptotic p-value may be non-finite, in which case the row is reported as `NA`.
**Note on Accessors**: EDI count-regression timings in this table explicitly call both `compute_wald_confidence_interval()` and `compute_wald_two_sided_pval()` so those rows remain Wald-only and cannot dispatch to bootstrap or likelihood-based fallback paths.
**Solver-Only Prebuilds**: Benchmark setup prebuilds exposed observed-data design matrices, reduced design matrices, strata IDs, and other fixed working inputs outside the timed region when the implementation exposes those hooks. The timed region then measures the full-inference kernel on those fixed inputs.
**Limitation**: Some canonical comparators only expose formula-based APIs rather than comparable low-level fit kernels. Those rows remain included, but their canonical timings may still contain formula/model-frame overhead beyond the numerical solver, variance, and p-value work itself.
**Timing Note**: All timings are medians over 10 warmed runs measured with adaptive batched `system.time`; paths below 0.01 ms use `microbenchmark(times = 2000)` instead.
**Timing P-Value**: `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.
**Row Highlighting**: Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons.

<table>
  <thead>
    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.16</td><td>stats</td><td>t.test(pool)</td><td>0.23</td><td>1.43x</td><td>0.0967</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.14</td><td>stats</td><td>wilcox.test</td><td>1.01</td><td>7.32x</td><td>1.06e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.43</td><td>stats</td><td>lm.fit(interact)+Wald</td><td>0.68</td><td>1.57x</td><td>1.98e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.06</td><td>stats</td><td>lm.fit+Wald</td><td>0.08</td><td>1.4x</td><td>5.18e-05</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>4.46</td><td>quantreg</td><td>rq+summary</td><td>3.79</td><td>0.85x</td><td>0.0301</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.47</td><td>MASS</td><td>rlm+summary</td><td>2.06</td><td>4.42x</td><td>1.96e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidCMH</td><td>incidence</td><td>0.16</td><td>stats</td><td>mantelhaen</td><td>2.03</td><td>12.85x</td><td>8.94e-10</td><td>***</td></tr>
    <tr><td>InferenceIncidExactFisher</td><td>incidence</td><td>1.90</td><td>stats</td><td>fisher.test</td><td>1.52</td><td>0.8x</td><td>0.00719</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.84</td><td>stats</td><td>glm+gcomp(RD)+Wald</td><td>2.98</td><td>3.57x</td><td>6.36e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.96</td><td>stats</td><td>glm+gcomp(RR)+Wald</td><td>2.92</td><td>3.05x</td><td>3.09e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>1.14</td><td>stats</td><td>glm.fit+Wald(log)</td><td>6.50</td><td>5.71x</td><td>3.89e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.32</td><td>stats</td><td>glm.fit+Wald</td><td>0.68</td><td>2.14x</td><td>5.2e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.18</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>1.15</td><td>6.55x</td><td>2.36e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.21</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>0.94</td><td>4.53x</td><td>2.36e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.97</td><td>stats</td><td>glm.fit(probit)+Wald</td><td>1.29</td><td>1.32x</td><td>0.00651</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.29</td><td>stats</td><td>prop.test</td><td>0.88</td><td>2.98x</td><td>0.000464</td><td>***</td></tr>
    <tr><td>InferenceCountHurdleNegBin</td><td>count</td><td>21.95</td><td>pscl</td><td>hurdle(nb)+summary</td><td>20.00</td><td>0.91x</td><td>0.347</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>2.98</td><td>pscl</td><td>hurdle+summary</td><td>11.74</td><td>3.94x</td><td>4.89e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>3.56</td><td>MASS</td><td>glm.nb+summary</td><td>18.92</td><td>5.32x</td><td>7.87e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.67</td><td>stats</td><td>glm.fit+Wald</td><td>0.96</td><td>1.44x</td><td>4.72e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.91</td><td>stats</td><td>glm.fit+Wald(quasi)</td><td>1.20</td><td>1.33x</td><td>0.0458</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.89</td><td>sandwich</td><td>glm+vcovHC</td><td>6.47</td><td>7.28x</td><td>3.17e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>2.89</td><td>pscl</td><td>zeroinfl(nb)+summary</td><td>34.64</td><td>12x</td><td>4.15e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>2.51</td><td>pscl</td><td>zeroinfl+summary</td><td>21.68</td><td>8.63x</td><td>2.94e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>2.21</td><td>betareg</td><td>betareg.fit+Wald</td><td>26.72</td><td>12.12x</td><td>6.29e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.88</td><td>stats</td><td>glm(quasi)+gcomp+Wald</td><td>3.17</td><td>3.61x</td><td>4.28e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.25</td><td>survival</td><td>coxph.fit(breslow)+Wald</td><td>0.48</td><td>1.96x</td><td>7.22e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>8.32</td><td>survival</td><td>coxph(null)+KM weighted residual mean diff + survdiff(rho=1)</td><td>10.44</td><td>1.25x</td><td>0.0255</td><td>*</td></tr>
    <tr><td>InferenceSurvivalKMDiff</td><td>survival</td><td>6.14</td><td>survival</td><td>survfit(median)+CI</td><td>6.08</td><td>0.99x</td><td>0.979</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>1.19</td><td>survival</td><td>survdiff</td><td>2.80</td><td>2.35x</td><td>1.66e-08</td><td>***</td></tr>
    <tr><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.60</td><td>survival</td><td>coxph.fit(strat,breslow)+Wald</td><td>0.53</td><td>0.88x</td><td>0.389</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.53</td><td>survival</td><td>survreg+summary</td><td>3.98</td><td>7.49x</td><td>2.22e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.78</td><td>VGAM</td><td>vglm+summary</td><td>20.54</td><td>26.47x</td><td>3.33e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.77</td><td>VGAM</td><td>vglm+summary</td><td>24.00</td><td>31.04x</td><td>2.76e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>1.14</td><td>ordinal</td><td>clm+gcomp+Wald</td><td>9.17</td><td>8.06x</td><td>6.1e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>0.29</td><td>clinfun</td><td>jonckheere</td><td>0.81</td><td>2.76x</td><td>1.01e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.89</td><td>ordinal</td><td>clm+summary</td><td>6.25</td><td>7.04x</td><td>2.65e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>0.10</td><td>stats</td><td>mean(ridit)</td><td>0.25</td><td>2.49x</td><td>2.3e-06</td><td>***</td></tr>
  </tbody>
</table>

<style>
    body, .markdown-body, .container {
        max-width: 1200px !important;
        width: 100% !important;
        margin: 0 auto !important;
    }
</style>
