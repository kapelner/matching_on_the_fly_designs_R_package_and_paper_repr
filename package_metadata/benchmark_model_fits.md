# EDI Exhaustive C++ Model Fit Benchmarks

_Generated: 2026-07-30 09:03:31 JST_

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Compilation Context

These rows are read from build metadata compiled into the loaded `EDI` shared object via `edi_build_info_cpp()`.

**Compilation warning:** EDI model-fit timings are sensitive to the compiler flags used to build the loaded `EDI.so`. If EDI is compiled without the proper optimized flags, or with flags that are known to degrade these kernels such as problematic LTO builds, the benchmark can show substantial performance regressions that reflect the binary build rather than the modeling algorithms.

*   **EDI shared object:** `/home/kapelner/R/x86_64-pc-linux-gnu-library/4.7/EDI/libs/EDI.so`
*   **EDI shared object mtime:** `2026-07-30 00:36:41`
*   **Capture method:** `configure-generated header compiled into EDI.so`
*   **Build timestamp:** `2026-07-29 12:12:22 JST`
*   **Build host:** `LAPTOP-J2T9TGGB`
*   **R version at build:** `R Under development (unstable) (2026-04-23 r89955) -- "Unsuffered Consequences"`
*   **R `CXX20` at build:** `g++`
*   **R `CXX20STD` at build:** `-std=gnu++20`
*   **R `CXX20FLAGS` at build:** `-O3 -march=native -funroll-loops -fno-math-errno -UNDEBUG -Wall -pedantic -g -O0`
*   **R `SHLIB_OPENMP_CXXFLAGS` at build:** `unavailable`
*   **Build env at build:** `EDI_PORTABLE=0`, `EDI_DISABLE_VECTORIZATION=0`, `EDI_NATIVE_SPEED=1`, `EDI_NATIVE_LTO=0`
*   **Package `PKG_CPPFLAGS` at build:** ``
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
*   **Averaging:** All timings are medians over 30 cold estimate-only timing samples measured with adaptive batched `system.time`; paths below 0.01 ms use `microbenchmark(times = 5000)` instead.
*   **Timing P-Value:** `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.
*   **Row Highlighting:** Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons from a failed fit; light blue rows are estimators with no canonical R implementation at all (only EDI is timed, `Canonical Time`/`Speedup`/`Timing Pval` are `NA` by design).
*   **Constraints**: Most matched-pair/KK and highly custom paths are excluded as per user request; the exceptions are the four light-blue rows below, whose custom joint likelihood or estimator family (KK combined matched+reservoir, Weibull frailty, zero-one-inflated beta) has no canonical R implementation to compare against — see `package_metadata/python_bindings_package_spec.md` for the underlying baseline-gap analysis.

## Results

<table>
  <thead>
    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.12</td><td>stats</td><td>HL median pairwise diff</td><td>1.38</td><td>11.34x</td><td>2.55e-23</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.02</td><td>stats</td><td>lm.fit</td><td>0.13</td><td>7.73x</td><td>2.22e-34</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>1.22</td><td>quantreg</td><td>rq.fit</td><td>1.21</td><td>0.99x</td><td>0.465</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.09</td><td>MASS</td><td>rlm(MM)</td><td>43.13</td><td>473.16x</td><td>6.22e-42</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.13</td><td>stats</td><td>glm.fit(ident)</td><td>11.26</td><td>88.59x</td><td>5.12e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.15</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>1.37</td><td>8.91x</td><td>1.02e-21</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.15</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>1.37</td><td>9x</td><td>4.63e-24</td><td>***</td></tr>
    <tr style="background-color: #cfe2ff;"><td>InferenceIncidKKCondLogitPlusGLMMOneLik</td><td>incidence</td><td>10.84</td><td>None</td><td>no canonical R implementation</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>2.52</td><td>stats</td><td>glm.fit(log)</td><td>11.60</td><td>4.6x</td><td>4.41e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.14</td><td>stats</td><td>glm.fit</td><td>1.58</td><td>11.05x</td><td>2.65e-22</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.09</td><td>stats</td><td>glm.fit(modified)</td><td>1.62</td><td>17.19x</td><td>4.34e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.30</td><td>stats</td><td>glm.fit(probit)</td><td>1.77</td><td>5.82x</td><td>4.25e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.02</td><td>stats</td><td>lm.fit(LPM)</td><td>0.09</td><td>5.79x</td><td>6.75e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>1.60</td><td>pscl</td><td>hurdle(nb)</td><td>41.20</td><td>25.68x</td><td>1.51e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>1.38</td><td>pscl</td><td>hurdle</td><td>16.57</td><td>12.01x</td><td>1.13e-26</td><td>***</td></tr>
    <tr style="background-color: #cfe2ff;"><td>InferenceCountKKCondPoissonOneLik</td><td>count</td><td>0.05</td><td>None</td><td>no canonical R implementation</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>0.41</td><td>MASS</td><td>glm.nb</td><td>49.20</td><td>118.88x</td><td>2.47e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.11</td><td>stats</td><td>glm.fit</td><td>1.91</td><td>17.79x</td><td>3.11e-17</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.11</td><td>stats</td><td>glm.fit(quasi)</td><td>1.56</td><td>14.85x</td><td>6.31e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.10</td><td>stats</td><td>glm.fit</td><td>1.65</td><td>16.41x</td><td>2.5e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>0.98</td><td>pscl</td><td>zeroinfl(nb)</td><td>133.00</td><td>136.05x</td><td>8.41e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>3.31</td><td>pscl</td><td>zeroinfl</td><td>62.25</td><td>18.78x</td><td>1.19e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>0.92</td><td>betareg</td><td>betareg.fit</td><td>27.50</td><td>29.95x</td><td>2.35e-30</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.13</td><td>stats</td><td>glm.fit(quasi)</td><td>1.14</td><td>8.53x</td><td>3.65e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.18</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>1.51</td><td>8.57x</td><td>2.57e-21</td><td>***</td></tr>
    <tr style="background-color: #cfe2ff;"><td>InferencePropZeroOneInflatedBetaRegr</td><td>proportion</td><td>3.09</td><td>None</td><td>no canonical R implementation</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.31</td><td>survival</td><td>coxph.fit(breslow)</td><td>0.56</td><td>1.77x</td><td>7.43e-15</td><td>***</td></tr>
    <tr style="background-color: #cfe2ff;"><td>InferenceSurvivalKKWeibullFrailtyOneLik</td><td>survival</td><td>1.33</td><td>None</td><td>no canonical R implementation</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.01</td><td>survival</td><td>survfit(median)</td><td>3.01</td><td>219.35x</td><td>8.95e-39</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.01</td><td>survival</td><td>survdiff</td><td>1.67</td><td>124.9x</td><td>2.47e-36</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.01</td><td>survival</td><td>survfit(rmean)</td><td>2.29</td><td>162.32x</td><td>4.12e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.52</td><td>survival</td><td>coxph.fit(strat)</td><td>0.67</td><td>1.28x</td><td>2.05e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.05</td><td>survival</td><td>survreg</td><td>3.34</td><td>64.15x</td><td>2.65e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.42</td><td>VGAM</td><td>vglm(acat)</td><td>12.92</td><td>30.66x</td><td>1.21e-31</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>0.44</td><td>ordinal</td><td>clm(cauchit)</td><td>9.00</td><td>20.44x</td><td>3.38e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>0.40</td><td>ordinal</td><td>clm(cloglog)</td><td>6.66</td><td>16.57x</td><td>6.45e-37</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.19</td><td>VGAM</td><td>vglm(cratio)</td><td>10.89</td><td>56x</td><td>5.15e-38</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.46</td><td>ordinal</td><td>clm+gcomp</td><td>11.73</td><td>25.39x</td><td>4.84e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>0.37</td><td>ordinal</td><td>clm(probit)</td><td>5.78</td><td>15.84x</td><td>4.59e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.43</td><td>ordinal</td><td>clm</td><td>6.08</td><td>13.99x</td><td>1.16e-36</td><td>***</td></tr>
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
**Row Highlighting**: Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons from a failed fit; light blue rows are estimators with no canonical R implementation at all (only EDI is timed, `Canonical Time`/`Speedup`/`Timing Pval` are `NA` by design).

<table>
  <thead>
    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.03</td><td>stats</td><td>t.test(pool)</td><td>0.12</td><td>4.05x</td><td>6.35e-37</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.04</td><td>stats</td><td>wilcox.test</td><td>0.54</td><td>15.24x</td><td>4.09e-37</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.20</td><td>stats</td><td>lm.fit(interact)+Wald</td><td>0.43</td><td>2.1x</td><td>5.27e-31</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.02</td><td>stats</td><td>lm.fit+Wald</td><td>0.06</td><td>3.92x</td><td>2.37e-37</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>1.73</td><td>quantreg</td><td>rq+summary</td><td>1.73</td><td>1x</td><td>0.25</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.04</td><td>MASS</td><td>rlm+summary</td><td>1.12</td><td>25.05x</td><td>4.75e-41</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidExactFisher</td><td>incidence</td><td>0.67</td><td>stats</td><td>fisher.test</td><td>0.76</td><td>1.13x</td><td>2e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.06</td><td>stats</td><td>glm+gcomp(RD)+Wald</td><td>1.68</td><td>27.56x</td><td>4.91e-41</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.06</td><td>stats</td><td>glm+gcomp(RR)+Wald</td><td>1.68</td><td>27.46x</td><td>4.36e-33</td><td>***</td></tr>
    <tr style="background-color: #cfe2ff;"><td>InferenceIncidKKCondLogitPlusGLMMOneLik</td><td>incidence</td><td>3.84</td><td>None</td><td>no canonical R implementation</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>0.03</td><td>stats</td><td>glm.fit+Wald(log)</td><td>3.94</td><td>149.54x</td><td>3.32e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.05</td><td>stats</td><td>glm.fit+Wald</td><td>0.57</td><td>10.98x</td><td>8.81e-34</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.01</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>0.57</td><td>55.34x</td><td>8.28e-41</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.08</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>0.67</td><td>7.9x</td><td>1.1e-37</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.10</td><td>stats</td><td>glm.fit(probit)+Wald</td><td>0.87</td><td>9.15x</td><td>2.54e-37</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.02</td><td>stats</td><td>prop.test</td><td>0.37</td><td>22.96x</td><td>1.02e-69</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>0.25</td><td>pscl</td><td>hurdle(nb)+summary</td><td>11.38</td><td>45.64x</td><td>4.83e-42</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>0.36</td><td>pscl</td><td>hurdle+summary</td><td>8.17</td><td>22.83x</td><td>2.47e-27</td><td>***</td></tr>
    <tr style="background-color: #cfe2ff;"><td>InferenceCountKKCondPoissonOneLik</td><td>count</td><td>0.04</td><td>None</td><td>no canonical R implementation</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>0.12</td><td>MASS</td><td>glm.nb+summary</td><td>8.65</td><td>72.13x</td><td>4.08e-38</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.05</td><td>stats</td><td>glm.fit+Wald</td><td>0.74</td><td>16.53x</td><td>6.78e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.04</td><td>stats</td><td>glm.fit+Wald(quasi)</td><td>0.74</td><td>19.42x</td><td>8.27e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.08</td><td>sandwich</td><td>glm+vcovHC</td><td>3.29</td><td>39.69x</td><td>1.08e-41</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>0.88</td><td>pscl</td><td>zeroinfl(nb)+summary</td><td>81.67</td><td>93.28x</td><td>3.36e-40</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>1.50</td><td>pscl</td><td>zeroinfl+summary</td><td>26.00</td><td>17.33x</td><td>9.87e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>0.21</td><td>betareg</td><td>betareg+summary</td><td>12.79</td><td>61.29x</td><td>1.39e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.06</td><td>stats</td><td>glm(quasi)+gcomp+Wald</td><td>1.65</td><td>27.57x</td><td>2.32e-38</td><td>***</td></tr>
    <tr style="background-color: #cfe2ff;"><td>InferencePropZeroOneInflatedBetaRegr</td><td>proportion</td><td>0.50</td><td>None</td><td>no canonical R implementation</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.11</td><td>survival</td><td>coxph.fit(breslow)+Wald</td><td>0.36</td><td>3.41x</td><td>2.91e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>1.23</td><td>survival</td><td>survdiff(rho=1)</td><td>1.30</td><td>1.06x</td><td>0.00117</td><td>**</td></tr>
    <tr style="background-color: #cfe2ff;"><td>InferenceSurvivalKKWeibullFrailtyOneLik</td><td>survival</td><td>0.72</td><td>None</td><td>no canonical R implementation</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>2.44</td><td>survival</td><td>survfit(median)+CI</td><td>2.63</td><td>1.08x</td><td>1.56e-06</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>9.80e-03</td><td>survival</td><td>survdiff</td><td>1.29</td><td>131.39x</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.38</td><td>survival</td><td>coxph.fit(strat)+Wald</td><td>0.51</td><td>1.37x</td><td>4.59e-31</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.06</td><td>survival</td><td>survreg+summary</td><td>2.54</td><td>45.12x</td><td>7.47e-37</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.13</td><td>VGAM</td><td>vglm+summary</td><td>11.53</td><td>86.05x</td><td>9.98e-35</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.06</td><td>VGAM</td><td>vglm+summary</td><td>13.17</td><td>208.73x</td><td>4.5e-16</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.32</td><td>ordinal</td><td>clm+gcomp+Wald</td><td>6.78</td><td>21.47x</td><td>5.87e-28</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>8.90e-03</td><td>clinfun</td><td>jonckheere</td><td>0.35</td><td>39.77x</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.14</td><td>ordinal</td><td>clm+summary</td><td>4.55</td><td>31.65x</td><td>5.89e-38</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>8.90e-03</td><td>stats</td><td>mean(ridit)</td><td>0.14</td><td>15.38x</td><td>NA</td><td></td></tr>
  </tbody>
</table>

## Utility / Math Kernel Performance

This table is an exhaustive inventory of every fast_* scalar math kernel in EDI/src — used inside the NegBin/Beta/ZINB/Hurdle likelihoods, KK21 negative-binomial fitting, probit cold-start heuristics, ordinal cloglog/cauchit link derivatives, and LRT/score-test p-values — against base R's vectorized equivalents.
Unlike the model-fit tables above, these are not full estimators: each row is a single vectorized special-function evaluation (digamma, trigamma, log-gamma, log-beta, error function, normal CDF/PDF/quantile (both closed-form and log-scale), chi-squared upper-tail p-value, arctangent, softplus, and the mu-parameterized negative-binomial density) over a fixed-length vector, with no design matrix, optimizer, or R6 object involved on either side.
**Vector length**: All rows evaluate the function over a vector of length $N=5000$; inputs are drawn from a domain realistic for each kernel's actual call sites in EDI (e.g. shape/rate-like values in $(0.5, 50)$ for digamma/trigamma/lgamma, probabilities in $(0, 1)$ for qnorm).
**Bare Metal EDI Timing**: Eight kernels (digamma, trigamma, lgamma, lbeta, dnbinom_mu, qnorm, log_pnorm, log_dnorm) call the package's own exported `fast_*_vec_cpp` wrapper (`EDI/src/fast_math_utils.cpp`, documented and exported via roxygen/NAMESPACE like any other EDI function) directly from the installed package — no separate compile step. The remaining six (pchisq_upper, erfc, pnorm_fast, dnorm_fast, atan, log1pexp) aren't promoted to package exports yet, so those rows call a thin bare-metal wrapper compiled standalone via `Rcpp::sourceCpp()` against the same `EDI/src` leaf headers (`benchmark/fast_math_utils_bench.cpp`) — same convention as `benchmark/fast_trigamma_speed_compare.cpp`. Either way, the wrapper just loops the internal `fast_*` scalar kernel over the input vector — no R6, no caching, no warm starts.
**Canonical R Timing**: The base R column calls the corresponding vectorized base/stats function directly on the same input vector (e.g. `digamma()`, `lgamma()`, `stats::qnorm()`, `stats::dnbinom(..., mu = , log = TRUE)`); `erfc` and `log1pexp` have no dedicated base R function, so those two rows use the standard base-R composition (`2*pnorm(-x*sqrt(2))`, `log1p(exp(x))`) a base-R user would otherwise write by hand.
**Timing Note**: All timings use the same adaptive batched `system.time`/`microbenchmark` harness as the tables above (medians over 30 cold samples; paths below 0.01 ms fall back to `microbenchmark(times = 5000)`).
**Timing P-Value**: `Timing Pval` reports a Welch two-sample t-test comparing the EDI and base R timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.
**Row Highlighting**: Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons from a failed evaluation.

<table>
  <thead>
    <tr><th>Function</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr style="background-color: #d9fdd3;"><td>dnorm_fast</td><td>utility</td><td>0.07</td><td>base/stats</td><td>dnorm</td><td>0.14</td><td>2x</td><td>4.08e-21</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_atan</td><td>utility</td><td>0.06</td><td>base/stats</td><td>atan</td><td>0.09</td><td>1.63x</td><td>4.36e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_digamma</td><td>utility</td><td>0.09</td><td>base/stats</td><td>digamma</td><td>0.47</td><td>5.45x</td><td>1.18e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_dnbinom_mu</td><td>utility</td><td>0.49</td><td>base/stats</td><td>dnbinom(mu=, log=TRUE)</td><td>0.66</td><td>1.35x</td><td>6.8e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_erfc</td><td>utility</td><td>0.12</td><td>base/stats</td><td>2*pnorm(-x*sqrt(2))</td><td>0.35</td><td>2.96x</td><td>6.5e-21</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_lbeta</td><td>utility</td><td>0.26</td><td>base/stats</td><td>lbeta</td><td>0.60</td><td>2.3x</td><td>8.72e-38</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_lgamma</td><td>utility</td><td>0.09</td><td>base/stats</td><td>lgamma</td><td>0.19</td><td>2.05x</td><td>5.21e-35</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_log1pexp</td><td>utility</td><td>0.10</td><td>base/stats</td><td>log1p(exp(x))</td><td>0.14</td><td>1.42x</td><td>1.59e-16</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_log_dnorm</td><td>utility</td><td>0.03</td><td>base/stats</td><td>dnorm(log=TRUE)</td><td>0.11</td><td>3.35x</td><td>1.34e-23</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_log_pnorm</td><td>utility</td><td>0.13</td><td>base/stats</td><td>pnorm(log.p=TRUE)</td><td>0.33</td><td>2.44x</td><td>1.13e-22</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_pchisq_upper</td><td>utility</td><td>0.82</td><td>base/stats</td><td>pchisq(lower.tail=FALSE)</td><td>1.05</td><td>1.28x</td><td>0.0442</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_qnorm</td><td>utility</td><td>0.06</td><td>base/stats</td><td>qnorm</td><td>0.13</td><td>2.23x</td><td>1.19e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>fast_trigamma</td><td>utility</td><td>0.05</td><td>base/stats</td><td>trigamma</td><td>0.65</td><td>13.41x</td><td>2.59e-31</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>pnorm_fast</td><td>utility</td><td>0.10</td><td>base/stats</td><td>pnorm</td><td>0.31</td><td>3.11x</td><td>1.71e-28</td><td>***</td></tr>
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
