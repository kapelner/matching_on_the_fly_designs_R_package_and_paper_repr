# EDI Exhaustive C++ Model Fit Benchmarks

_Generated: 2026-07-01 20:41:11 JST_

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Compilation Context

These rows are read from build metadata compiled into the loaded `EDI` shared object via `edi_build_info_cpp()`.

**Compilation warning:** EDI model-fit timings are sensitive to the compiler flags used to build the loaded `EDI.so`. If EDI is compiled without the proper optimized flags, or with flags that are known to degrade these kernels such as problematic LTO builds, the benchmark can show substantial performance regressions that reflect the binary build rather than the modeling algorithms.

*   **EDI shared object:** `/home/kapelner/R/x86_64-pc-linux-gnu-library/4.7/EDI/libs/EDI.so`
*   **EDI shared object mtime:** `2026-07-01 09:33:33`
*   **Capture method:** `configure-generated header compiled into EDI.so`
*   **Build timestamp:** `2026-06-29 14:02:11 JST`
*   **Build host:** `LAPTOP-J2T9TGGB`
*   **R version at build:** `R Under development (unstable) (2026-04-23 r89955) -- "Unsuffered Consequences"`
*   **R `CXX20` at build:** `g++`
*   **R `CXX20STD` at build:** `-std=gnu++20`
*   **R `CXX20FLAGS` at build:** `-O3 -march=native -funroll-loops -fno-math-errno -UNDEBUG -Wall -pedantic -g -O0`
*   **R `SHLIB_OPENMP_CXXFLAGS` at build:** `unavailable`
*   **Build env at build:** `EDI_PORTABLE=0`, `EDI_DISABLE_VECTORIZATION=0`, `EDI_NATIVE_SPEED=1`, `EDI_NATIVE_LTO=0`
*   **Package `PKG_CPPFLAGS` at build:** `-I../inst/include`
*   **Package `PKG_CXXFLAGS` at build:** `$(SHLIB_OPENMP_CXXFLAGS) -DNDEBUG -DEIGEN_NO_DEBUG -Wno-ignored-attributes -march=native -mtune=native -fno-lto; override CXXFLAGS+=-O3`
*   **Package `PKG_LIBS` at build:** `$(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -ltbb12 -fstack-protector`
*   **Compiler reported by binary:** `15.2.0`
*   **Compiler optimization macro enabled:** `TRUE`
*   **Compiler fast-math macro enabled:** `FALSE`
*   **Eigen vectorization disabled macro enabled:** `TRUE`

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
*   **Averaging:** All timings are medians over 30 cold estimate-only timing samples measured with adaptive batched `system.time`; paths below 0.01 ms use `microbenchmark(times = 5000)` instead.
*   **Timing P-Value:** `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.
*   **Row Highlighting:** Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons.
*   **Constraints**: Matched-pair/KK and highly custom paths are excluded as per user request.

## Results

<table>
  <thead>
    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.60</td><td>stats</td><td>HL median pairwise diff</td><td>1.32</td><td>2.19x</td><td>1.24e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.03</td><td>stats</td><td>lm.fit</td><td>0.09</td><td>3.36x</td><td>2.13e-21</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>1.46</td><td>quantreg</td><td>rq.fit</td><td>1.43</td><td>0.98x</td><td>0.0441</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.34</td><td>MASS</td><td>rlm(MM)</td><td>44.20</td><td>131.13x</td><td>4.87e-21</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.12</td><td>stats</td><td>glm.fit(ident)</td><td>10.91</td><td>90.26x</td><td>3.95e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.18</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>1.44</td><td>8.12x</td><td>3.91e-20</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.20</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>1.51</td><td>7.61x</td><td>5.52e-17</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>1.52</td><td>stats</td><td>glm.fit(log)</td><td>3.98</td><td>2.61x</td><td>7.91e-23</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.13</td><td>stats</td><td>glm.fit</td><td>1.51</td><td>11.68x</td><td>1.8e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.24</td><td>stats</td><td>glm.fit(modified)</td><td>1.68</td><td>7.1x</td><td>4.15e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.53</td><td>stats</td><td>glm.fit(probit)</td><td>1.60</td><td>3.04x</td><td>1.02e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.02</td><td>stats</td><td>lm.fit(LPM)</td><td>0.09</td><td>3.78x</td><td>2.65e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>1.73</td><td>pscl</td><td>hurdle(nb)</td><td>39.83</td><td>23.02x</td><td>3.16e-49</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>1.79</td><td>pscl</td><td>hurdle</td><td>20.71</td><td>11.57x</td><td>2.35e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>0.46</td><td>MASS</td><td>glm.nb</td><td>48.20</td><td>105.55x</td><td>2.64e-50</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.23</td><td>stats</td><td>glm.fit</td><td>1.65</td><td>7.33x</td><td>4.95e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.22</td><td>stats</td><td>glm.fit(quasi)</td><td>1.63</td><td>7.34x</td><td>8.78e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.27</td><td>stats</td><td>glm.fit</td><td>1.41</td><td>5.29x</td><td>9.98e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>1.73</td><td>pscl</td><td>zeroinfl(nb)</td><td>143.75</td><td>83.07x</td><td>2.41e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>3.95</td><td>pscl</td><td>zeroinfl</td><td>69.83</td><td>17.66x</td><td>6.53e-22</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>1.69</td><td>betareg</td><td>betareg.fit</td><td>30.50</td><td>18.1x</td><td>1.48e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.13</td><td>stats</td><td>glm.fit(quasi)</td><td>1.18</td><td>9.39x</td><td>1.36e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.15</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>1.27</td><td>8.7x</td><td>5.66e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.35</td><td>survival</td><td>coxph.fit(breslow)</td><td>0.55</td><td>1.55x</td><td>6.07e-18</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.01</td><td>survival</td><td>survfit(median)</td><td>3.23</td><td>224.57x</td><td>3.45e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.02</td><td>survival</td><td>survdiff</td><td>1.71</td><td>81.22x</td><td>8.24e-45</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.01</td><td>survival</td><td>survfit(rmean)</td><td>2.60</td><td>181.46x</td><td>8.18e-20</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.59</td><td>survival</td><td>coxph.fit(strat)</td><td>0.69</td><td>1.17x</td><td>3e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.10</td><td>survival</td><td>survreg</td><td>3.72</td><td>39.15x</td><td>8.34e-30</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.52</td><td>VGAM</td><td>vglm(acat)</td><td>12.11</td><td>23.11x</td><td>9.01e-30</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>0.70</td><td>ordinal</td><td>clm(cauchit)</td><td>8.07</td><td>11.59x</td><td>4.01e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>0.62</td><td>ordinal</td><td>clm(cloglog)</td><td>6.80</td><td>10.93x</td><td>5.03e-38</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.43</td><td>VGAM</td><td>vglm(cratio)</td><td>11.32</td><td>26.04x</td><td>2.16e-34</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.92</td><td>ordinal</td><td>clm+gcomp</td><td>12.09</td><td>13.2x</td><td>8.37e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>0.72</td><td>ordinal</td><td>clm(probit)</td><td>6.16</td><td>8.58x</td><td>1.08e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.87</td><td>ordinal</td><td>clm</td><td>7.53</td><td>8.67x</td><td>4.28e-22</td><td>***</td></tr>
  </tbody>
</table>

## Wald Test Performance (Full Inference)

This table compares the performance of **Full Inference** (Model Fit + Standard Error calculation + P-value derivation).
Unlike the point-estimation table above, these results include the computational cost of the variance-covariance matrix (Hessian or Fisher Information) and the Wald test statistic calculation.
All paths (EDI and Canonical) use a reduced sample size ($N=200$) for this full-inference benchmark to ensure iterative stability.
**Stratified Cox Exception**: For `InferenceSurvivalStratCoxPHRegr`, the benchmark injects low-cardinality covariates before outcome generation so the row exercises a genuinely stratified Cox fit rather than the unstratified fallback.
EDI regression models (Logistic, Poisson) are benchmarked using the **IRLS** optimizer for these Wald tests.
**Solver-Only Prebuilds**: Benchmark setup prebuilds exposed observed-data design matrices, reduced design matrices, strata IDs, and other fixed working inputs outside the timed region when the implementation exposes those hooks. The timed region then measures the full-inference kernel on those fixed inputs.
**Limitation**: Some canonical comparators only expose formula-based APIs rather than comparable low-level fit kernels. Those rows remain included, but their canonical timings may still contain formula/model-frame overhead beyond the numerical solver, variance, and p-value work itself.
**Timing Note**: All timings are medians over 30 warmed runs measured with adaptive batched `system.time`; paths below 0.01 ms use `microbenchmark(times = 5000)` instead.
**Timing P-Value**: `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.
**Row Highlighting**: Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons.

<table>
  <thead>
    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.03</td><td>stats</td><td>t.test(pool)</td><td>0.12</td><td>4.22x</td><td>7.45e-41</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.08</td><td>stats</td><td>wilcox.test</td><td>0.66</td><td>8.28x</td><td>6.99e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.18</td><td>stats</td><td>lm.fit(interact)+Wald</td><td>0.42</td><td>2.34x</td><td>2.77e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.01</td><td>stats</td><td>lm.fit+Wald</td><td>0.06</td><td>4.3x</td><td>2.03e-24</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>2.49</td><td>quantreg</td><td>rq+summary</td><td>1.86</td><td>0.75x</td><td>2.49e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.08</td><td>MASS</td><td>rlm+summary</td><td>1.10</td><td>13.52x</td><td>3.38e-36</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidExactFisher</td><td>incidence</td><td>0.65</td><td>stats</td><td>fisher.test</td><td>0.73</td><td>1.13x</td><td>4.47e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.07</td><td>stats</td><td>glm+gcomp(RD)+Wald</td><td>1.80</td><td>24.15x</td><td>9.25e-30</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.07</td><td>stats</td><td>glm+gcomp(RR)+Wald</td><td>1.69</td><td>23.19x</td><td>2.2e-22</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>1.13</td><td>stats</td><td>glm.fit+Wald(log)</td><td>2.97</td><td>2.63x</td><td>5.99e-42</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.04</td><td>stats</td><td>glm.fit+Wald</td><td>0.59</td><td>13.49x</td><td>7.89e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.01</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>0.57</td><td>55.85x</td><td>4.35e-39</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.08</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>0.68</td><td>8.61x</td><td>5.22e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.14</td><td>stats</td><td>glm.fit(probit)+Wald</td><td>0.75</td><td>5.37x</td><td>9.34e-71</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.01</td><td>stats</td><td>prop.test</td><td>0.35</td><td>25.57x</td><td>4.37e-45</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>0.26</td><td>pscl</td><td>hurdle(nb)+summary</td><td>10.43</td><td>40.63x</td><td>8.21e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>0.22</td><td>pscl</td><td>hurdle+summary</td><td>8.52</td><td>38.3x</td><td>2.1e-21</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>0.28</td><td>MASS</td><td>glm.nb+summary</td><td>11.81</td><td>42.44x</td><td>1.1e-39</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.07</td><td>stats</td><td>glm.fit+Wald</td><td>0.80</td><td>10.95x</td><td>1.8e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.08</td><td>stats</td><td>glm.fit+Wald(quasi)</td><td>0.87</td><td>11.37x</td><td>4.14e-34</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.12</td><td>sandwich</td><td>glm+vcovHC</td><td>3.46</td><td>29.32x</td><td>8.32e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>1.46</td><td>pscl</td><td>zeroinfl(nb)+summary</td><td>229.50</td><td>156.9x</td><td>1.61e-44</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>2.47</td><td>pscl</td><td>zeroinfl+summary</td><td>35.21</td><td>14.25x</td><td>1.51e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>0.45</td><td>betareg</td><td>betareg+summary</td><td>12.94</td><td>28.51x</td><td>6.78e-37</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.07</td><td>stats</td><td>glm(quasi)+gcomp+Wald</td><td>1.68</td><td>24.34x</td><td>1.53e-40</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.12</td><td>survival</td><td>coxph.fit(breslow)+Wald</td><td>0.40</td><td>3.23x</td><td>3.66e-48</td><td>***</td></tr>
    <tr><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>1.32</td><td>survival</td><td>survdiff(rho=1)</td><td>1.33</td><td>1.01x</td><td>0.858</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>2.32</td><td>survival</td><td>survfit(median)+CI</td><td>2.54</td><td>1.1x</td><td>1.99e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.01</td><td>survival</td><td>survdiff</td><td>1.29</td><td>123.91x</td><td>1.62e-45</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.39</td><td>survival</td><td>coxph.fit(strat)+Wald</td><td>0.60</td><td>1.53x</td><td>2.37e-15</td><td>***</td></tr>
    <tr><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.16</td><td>survival</td><td>survreg+summary</td><td>3.31</td><td>20.69x</td><td>0.33</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.15</td><td>VGAM</td><td>vglm+summary</td><td>13.33</td><td>87.25x</td><td>1.16e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.11</td><td>VGAM</td><td>vglm+summary</td><td>11.17</td><td>101.61x</td><td>1.1e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.43</td><td>ordinal</td><td>clm+gcomp+Wald</td><td>5.96</td><td>13.75x</td><td>1.59e-19</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>0.03</td><td>clinfun</td><td>jonckheere</td><td>0.36</td><td>12.4x</td><td>3.5e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.24</td><td>ordinal</td><td>clm+summary</td><td>4.85</td><td>20.08x</td><td>3.6e-36</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>0.01</td><td>stats</td><td>mean(ridit)</td><td>0.13</td><td>10.84x</td><td>3.26e-39</td><td>***</td></tr>
  </tbody>
</table>

## Garbage Collection and Cache Management

To ensure that the benchmark results are highly precise, reproducible, and represent the actual computation speed of the numerical solvers, the benchmarking harness uses the following garbage collection and cache management strategies:

### 1. Garbage Collection (GC) Filtering
Garbage collection cycles run automatically by the R interpreter and can introduce significant, arbitrary pauses that skew timing measurements. To isolate the execution time of the code from R's GC overhead:
* **GC Disabling**: We disable R's memory stress-testing mode using `gctorture(FALSE)` before running timing loops.
* **Proactive Compaction**: In the `system.time()` path, we invoke `gc(verbose = FALSE)` immediately before timing each replicate. This starts the timer on a clean, compacted heap, minimizing the likelihood of triggering an automatic garbage collection cycle mid-replicate.
* **Automatic Filtering**: In the microbenchmarking path, we utilize the `bench::mark()` engine with the `filter_gc = TRUE` parameter, which automatically tracks and discards timing iterations during which a garbage collection event occurred.

### 2. Cold-Start Guarantee for EDI and Symmetric Warm-Up for Both Sides
Both EDI and canonical timing expressions receive a single **validation/warm-up call** executed once before the calibration loop begins. This puts the machine code and working data into the instruction and data caches in the same warmed state for both sides, so the official timed replicates start on equal footing.

EDI timings call exported C++ functions directly — no R6 objects are instantiated during benchmarking. As a result, **no R6 result caches exist to manage**. Each call to the C++ solver (e.g. `fast_logistic_regression_cpp`, `fast_ordinal_regression_cpp`) starts from a freshly zero-initialized parameter vector (or a model-specific data-driven initialization when `smart_cold_start = TRUE`). No prior-fit results are carried across timing repetitions, so every replication is a genuine cold start for the numerical optimizer.

<style>
    body, .markdown-body, .container {
        max-width: 1200px !important;
        width: 100% !important;
        margin: 0 auto !important;
    }
</style>
