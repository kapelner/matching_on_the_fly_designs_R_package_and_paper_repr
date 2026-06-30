# EDI Exhaustive C++ Model Fit Benchmarks

_Generated: 2026-07-01 06:37:27 JST_

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Compilation Context

These rows are read from build metadata compiled into the loaded `EDI` shared object via `edi_build_info_cpp()`.

**Compilation warning:** EDI model-fit timings are sensitive to the compiler flags used to build the loaded `EDI.so`. If EDI is compiled without the proper optimized flags, or with flags that are known to degrade these kernels such as problematic LTO builds, the benchmark can show substantial performance regressions that reflect the binary build rather than the modeling algorithms.

*   **EDI shared object:** `/home/kapelner/R/x86_64-pc-linux-gnu-library/4.7/EDI/libs/EDI.so`
*   **EDI shared object mtime:** `2026-07-01 06:36:16`
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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.60</td><td>stats</td><td>HL median pairwise diff</td><td>1.30</td><td>2.16x</td><td>4.48e-17</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.03</td><td>stats</td><td>lm.fit</td><td>0.09</td><td>3.35x</td><td>7.87e-12</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>1.49</td><td>quantreg</td><td>rq.fit</td><td>1.49</td><td>1x</td><td>0.995</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.28</td><td>MASS</td><td>rlm(MM)</td><td>42.80</td><td>153.58x</td><td>2.95e-30</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.15</td><td>stats</td><td>glm.fit(ident)</td><td>11.70</td><td>79.29x</td><td>1.01e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.17</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>1.46</td><td>8.51x</td><td>4.64e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.15</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>1.38</td><td>9.22x</td><td>3.62e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>2.16</td><td>stats</td><td>glm.fit(log)</td><td>4.19</td><td>1.94x</td><td>2.21e-19</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.13</td><td>stats</td><td>glm.fit</td><td>1.42</td><td>10.94x</td><td>8.97e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.24</td><td>stats</td><td>glm.fit(modified)</td><td>1.77</td><td>7.27x</td><td>1.1e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.56</td><td>stats</td><td>glm.fit(probit)</td><td>1.68</td><td>3x</td><td>1.2e-23</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.03</td><td>stats</td><td>lm.fit(LPM)</td><td>0.09</td><td>3.51x</td><td>7.44e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>1.72</td><td>pscl</td><td>hurdle(nb)</td><td>39.58</td><td>23.02x</td><td>5.05e-41</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>1.59</td><td>pscl</td><td>hurdle</td><td>17.46</td><td>10.96x</td><td>6.02e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>0.46</td><td>MASS</td><td>glm.nb</td><td>47.00</td><td>103.18x</td><td>2.49e-42</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.22</td><td>stats</td><td>glm.fit</td><td>1.63</td><td>7.53x</td><td>2.07e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.22</td><td>stats</td><td>glm.fit(quasi)</td><td>1.60</td><td>7.23x</td><td>3.28e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.27</td><td>stats</td><td>glm.fit</td><td>1.47</td><td>5.45x</td><td>1.23e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>1.69</td><td>pscl</td><td>zeroinfl(nb)</td><td>138.50</td><td>81.89x</td><td>2.66e-31</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>4.01</td><td>pscl</td><td>zeroinfl</td><td>63.25</td><td>15.78x</td><td>1.37e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>1.56</td><td>betareg</td><td>betareg.fit</td><td>28.12</td><td>17.98x</td><td>8.07e-30</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.13</td><td>stats</td><td>glm.fit(quasi)</td><td>1.24</td><td>9.6x</td><td>2.07e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.15</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>1.33</td><td>8.69x</td><td>8.02e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.34</td><td>survival</td><td>coxph.fit(breslow)</td><td>0.53</td><td>1.57x</td><td>8.25e-41</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.01</td><td>survival</td><td>survfit(median)</td><td>3.34</td><td>233.56x</td><td>2.07e-43</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.02</td><td>survival</td><td>survdiff</td><td>1.90</td><td>86.46x</td><td>1.62e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.01</td><td>survival</td><td>survfit(rmean)</td><td>2.26</td><td>155.84x</td><td>2.08e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.61</td><td>survival</td><td>coxph.fit(strat)</td><td>0.70</td><td>1.15x</td><td>6.67e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.10</td><td>survival</td><td>survreg</td><td>3.42</td><td>35.57x</td><td>3.78e-38</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.75</td><td>VGAM</td><td>vglm(acat)</td><td>13.22</td><td>17.56x</td><td>3.76e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>0.71</td><td>ordinal</td><td>clm(cauchit)</td><td>8.54</td><td>12.01x</td><td>2.2e-41</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>0.64</td><td>ordinal</td><td>clm(cloglog)</td><td>7.69</td><td>11.92x</td><td>6.81e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.47</td><td>VGAM</td><td>vglm(cratio)</td><td>12.03</td><td>25.81x</td><td>3.59e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.70</td><td>ordinal</td><td>clm+gcomp</td><td>12.17</td><td>17.31x</td><td>1.77e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>0.75</td><td>ordinal</td><td>clm(probit)</td><td>6.65</td><td>8.81x</td><td>5.61e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.58</td><td>ordinal</td><td>clm</td><td>6.28</td><td>10.88x</td><td>1.27e-34</td><td>***</td></tr>
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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.03</td><td>stats</td><td>t.test(pool)</td><td>0.12</td><td>3.93x</td><td>5.93e-54</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.08</td><td>stats</td><td>wilcox.test</td><td>0.54</td><td>6.94x</td><td>8.78e-37</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.18</td><td>stats</td><td>lm.fit(interact)+Wald</td><td>0.45</td><td>2.53x</td><td>3.18e-49</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.01</td><td>stats</td><td>lm.fit+Wald</td><td>0.07</td><td>4.54x</td><td>2.2e-32</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>1.86</td><td>quantreg</td><td>rq+summary</td><td>1.87</td><td>1.01x</td><td>0.922</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.08</td><td>MASS</td><td>rlm+summary</td><td>1.25</td><td>14.71x</td><td>4.32e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidExactFisher</td><td>incidence</td><td>0.68</td><td>stats</td><td>fisher.test</td><td>0.75</td><td>1.1x</td><td>9.14e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.07</td><td>stats</td><td>glm+gcomp(RD)+Wald</td><td>1.85</td><td>25.19x</td><td>2.82e-38</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.07</td><td>stats</td><td>glm+gcomp(RR)+Wald</td><td>1.71</td><td>22.98x</td><td>4.94e-40</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>1.50</td><td>stats</td><td>glm.fit+Wald(log)</td><td>3.01</td><td>2x</td><td>2.46e-20</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.05</td><td>stats</td><td>glm.fit+Wald</td><td>0.61</td><td>13.33x</td><td>4.62e-38</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.01</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>0.59</td><td>57.19x</td><td>9.15e-41</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.09</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>0.71</td><td>7.87x</td><td>1.26e-30</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.14</td><td>stats</td><td>glm.fit(probit)+Wald</td><td>0.79</td><td>5.68x</td><td>2.7e-50</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.01</td><td>stats</td><td>prop.test</td><td>0.36</td><td>25.24x</td><td>2.45e-44</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>0.26</td><td>pscl</td><td>hurdle(nb)+summary</td><td>9.90</td><td>38.05x</td><td>2.21e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>0.23</td><td>pscl</td><td>hurdle+summary</td><td>7.73</td><td>33.6x</td><td>1.26e-44</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>0.28</td><td>MASS</td><td>glm.nb+summary</td><td>13.55</td><td>49.26x</td><td>4.13e-31</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.08</td><td>stats</td><td>glm.fit+Wald</td><td>0.88</td><td>11.12x</td><td>6.72e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.08</td><td>stats</td><td>glm.fit+Wald(quasi)</td><td>1.01</td><td>11.96x</td><td>1.06e-31</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.13</td><td>sandwich</td><td>glm+vcovHC</td><td>5.29</td><td>41.78x</td><td>8.13e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>1.59</td><td>pscl</td><td>zeroinfl(nb)+summary</td><td>251.00</td><td>157.87x</td><td>3.67e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>1.96</td><td>pscl</td><td>zeroinfl+summary</td><td>39.30</td><td>20.08x</td><td>1.57e-21</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>0.46</td><td>betareg</td><td>betareg+summary</td><td>14.70</td><td>31.61x</td><td>9.58e-30</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.07</td><td>stats</td><td>glm(quasi)+gcomp+Wald</td><td>1.81</td><td>25.04x</td><td>1.8e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.12</td><td>survival</td><td>coxph.fit(breslow)+Wald</td><td>0.38</td><td>3.02x</td><td>3.12e-40</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>1.32</td><td>survival</td><td>survdiff(rho=1)</td><td>1.40</td><td>1.06x</td><td>1.42e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>2.45</td><td>survival</td><td>survfit(median)+CI</td><td>2.66</td><td>1.09x</td><td>6.15e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.01</td><td>survival</td><td>survdiff</td><td>1.33</td><td>124.29x</td><td>6.08e-39</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.41</td><td>survival</td><td>coxph.fit(strat)+Wald</td><td>0.56</td><td>1.39x</td><td>3.4e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.09</td><td>survival</td><td>survreg+summary</td><td>2.77</td><td>31.83x</td><td>2.78e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.19</td><td>VGAM</td><td>vglm+summary</td><td>11.81</td><td>62.24x</td><td>4.21e-34</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.10</td><td>VGAM</td><td>vglm+summary</td><td>10.72</td><td>103.11x</td><td>9.6e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.44</td><td>ordinal</td><td>clm+gcomp+Wald</td><td>5.94</td><td>13.53x</td><td>5.83e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>0.03</td><td>clinfun</td><td>jonckheere</td><td>0.40</td><td>12.4x</td><td>3.88e-31</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.26</td><td>ordinal</td><td>clm+summary</td><td>4.82</td><td>18.69x</td><td>2.65e-38</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>0.01</td><td>stats</td><td>mean(ridit)</td><td>0.12</td><td>10.69x</td><td>1.08e-42</td><td>***</td></tr>
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
