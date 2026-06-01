# EDI Exhaustive C++ Model Fit Benchmarks

_Generated: 2026-06-01 02:04:29 EDT_

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Compilation Context

These rows are read from build metadata compiled into the loaded `EDI` shared object via `edi_build_info_cpp()`.

**Compilation warning:** EDI model-fit timings are sensitive to the compiler flags used to build the loaded `EDI.so`. If EDI is compiled without the proper optimized flags, or with flags that are known to degrade these kernels such as problematic LTO builds, the benchmark can show substantial performance regressions that reflect the binary build rather than the modeling algorithms.

*   **EDI shared object:** `/home/kapelner/R/x86_64-pc-linux-gnu-library/4.7/EDI/libs/EDI.so`
*   **EDI shared object mtime:** `2026-06-01 01:40:42`
*   **Capture method:** `configure-generated header compiled into EDI.so`
*   **Build timestamp:** `2026-05-31 19:49:31 EDT`
*   **Build host:** `LAPTOP-J2T9TGGB`
*   **R version at build:** `R Under development (unstable) (2026-04-23 r89955) -- "Unsuffered Consequences"`
*   **R `CXX20` at build:** `g++ -O3 -march=native -flto -fno-math-errno`
*   **R `CXX20STD` at build:** `-std=gnu++20`
*   **R `CXX20FLAGS` at build:** `-g -O2 -UNDEBUG -Wall -pedantic -g -O0`
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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.60</td><td>stats</td><td>HL median pairwise diff</td><td>2.01</td><td>3.37x</td><td>1.82e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.03</td><td>stats</td><td>lm.fit</td><td>0.14</td><td>4.92x</td><td>5.28e-10</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>1.74</td><td>quantreg</td><td>rq.fit</td><td>1.74</td><td>1x</td><td>0.787</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.30</td><td>MASS</td><td>rlm(MM)</td><td>42.90</td><td>145.3x</td><td>3.34e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.17</td><td>stats</td><td>glm.fit(ident)</td><td>14.78</td><td>84.57x</td><td>1.2e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.20</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>2.23</td><td>10.94x</td><td>7.61e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.19</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>1.92</td><td>10.06x</td><td>4.94e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>9.19</td><td>stats</td><td>glm.fit(log)</td><td>12.47</td><td>1.36x</td><td>0.00036</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.17</td><td>stats</td><td>glm.fit</td><td>2.35</td><td>14x</td><td>2.41e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.21</td><td>stats</td><td>glm.fit(modified)</td><td>2.08</td><td>9.71x</td><td>5.22e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.70</td><td>stats</td><td>glm.fit(probit)</td><td>2.07</td><td>2.96x</td><td>5.99e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.03</td><td>stats</td><td>lm.fit(LPM)</td><td>0.14</td><td>5.24x</td><td>1.69e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>8.72</td><td>pscl</td><td>hurdle(nb)</td><td>48.10</td><td>5.52x</td><td>1.08e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>1.93</td><td>pscl</td><td>hurdle</td><td>20.21</td><td>10.49x</td><td>1.41e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>5.93</td><td>MASS</td><td>glm.nb</td><td>49.20</td><td>8.29x</td><td>2.72e-18</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.24</td><td>stats</td><td>glm.fit</td><td>2.43</td><td>10.28x</td><td>4.06e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.24</td><td>stats</td><td>glm.fit(quasi)</td><td>2.31</td><td>9.42x</td><td>6.21e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.22</td><td>stats</td><td>glm.fit</td><td>1.96</td><td>8.76x</td><td>2.67e-16</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>2.34</td><td>pscl</td><td>zeroinfl(nb)</td><td>179.25</td><td>76.65x</td><td>9.72e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>5.41</td><td>pscl</td><td>zeroinfl</td><td>78.83</td><td>14.58x</td><td>3.39e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>5.82</td><td>betareg</td><td>betareg.fit</td><td>32.19</td><td>5.53x</td><td>3.75e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.17</td><td>stats</td><td>glm.fit(quasi)</td><td>1.70</td><td>9.93x</td><td>4.79e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.17</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>1.79</td><td>10.37x</td><td>4.2e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.41</td><td>survival</td><td>coxph.fit(breslow)</td><td>0.61</td><td>1.5x</td><td>2.45e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.02</td><td>survival</td><td>survfit(median)</td><td>4.31</td><td>284.5x</td><td>1.74e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.02</td><td>survival</td><td>survdiff</td><td>2.09</td><td>90.95x</td><td>1.33e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.02</td><td>survival</td><td>survfit(rmean)</td><td>2.64</td><td>167.11x</td><td>1.04e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.75</td><td>survival</td><td>coxph.fit(strat)</td><td>0.88</td><td>1.16x</td><td>0.000328</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.10</td><td>survival</td><td>survreg</td><td>4.10</td><td>41.24x</td><td>8.16e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.60</td><td>VGAM</td><td>vglm(acat)</td><td>16.20</td><td>26.96x</td><td>1.61e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>0.63</td><td>ordinal</td><td>clm(cauchit)</td><td>9.85</td><td>15.52x</td><td>5.02e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>0.93</td><td>ordinal</td><td>clm(cloglog)</td><td>9.20</td><td>9.89x</td><td>7.75e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.49</td><td>VGAM</td><td>vglm(cratio)</td><td>15.07</td><td>31.02x</td><td>1.4e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.85</td><td>ordinal</td><td>clm+gcomp</td><td>14.46</td><td>17.05x</td><td>1.79e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>0.80</td><td>ordinal</td><td>clm(probit)</td><td>8.64</td><td>10.78x</td><td>1.13e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.72</td><td>ordinal</td><td>clm</td><td>7.84</td><td>10.91x</td><td>5.4e-11</td><td>***</td></tr>
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


## Wald Test Performance (Full Inference)

This table compares the performance of **Full Inference** (Model Fit + Standard Error calculation + P-value derivation).
Unlike the point-estimation table above, these results include the computational cost of the variance-covariance matrix (Hessian or Fisher Information) and the Wald test statistic calculation.
All paths (EDI and Canonical) use a reduced sample size ($N=200$) for this full-inference benchmark to ensure iterative stability.
EDI timings in this table correspond to fixed `iBCRD` design objects.
**Stratified Cox Exception**: For `InferenceSurvivalStratCoxPHRegr`, the benchmark injects low-cardinality covariates before outcome generation so the row exercises a genuinely stratified Cox fit rather than the unstratified fallback.
EDI regression models (Logistic, Poisson) are benchmarked using the **IRLS** optimizer for these Wald tests.
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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.04</td><td>stats</td><td>t.test(pool)</td><td>0.16</td><td>3.77x</td><td>1.24e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.12</td><td>stats</td><td>wilcox.test</td><td>0.68</td><td>5.77x</td><td>4.2e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.22</td><td>stats</td><td>lm.fit(interact)+Wald</td><td>0.49</td><td>2.23x</td><td>7.61e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.02</td><td>stats</td><td>lm.fit+Wald</td><td>0.07</td><td>3.39x</td><td>8.06e-14</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>2.23</td><td>quantreg</td><td>rq+summary</td><td>2.16</td><td>0.97x</td><td>0.26</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.10</td><td>MASS</td><td>rlm+summary</td><td>1.25</td><td>12.53x</td><td>1.02e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidExactFisher</td><td>incidence</td><td>0.72</td><td>stats</td><td>fisher.test</td><td>0.88</td><td>1.23x</td><td>0.000124</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.10</td><td>stats</td><td>glm+gcomp(RD)+Wald</td><td>3.11</td><td>32.07x</td><td>2.78e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.10</td><td>stats</td><td>glm+gcomp(RR)+Wald</td><td>2.44</td><td>25.57x</td><td>8.2e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>0.80</td><td>stats</td><td>glm.fit+Wald(log)</td><td>5.82</td><td>7.28x</td><td>4.54e-17</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.07</td><td>stats</td><td>glm.fit+Wald</td><td>0.69</td><td>9.47x</td><td>1.41e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.01</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>0.71</td><td>57.63x</td><td>1.69e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.11</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>0.90</td><td>8.56x</td><td>1.87e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.18</td><td>stats</td><td>glm.fit(probit)+Wald</td><td>0.89</td><td>5.1x</td><td>1.37e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.02</td><td>stats</td><td>prop.test</td><td>0.40</td><td>19.42x</td><td>7.77e-16</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>1.72</td><td>pscl</td><td>hurdle(nb)+summary</td><td>16.75</td><td>9.74x</td><td>3.77e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>1.32</td><td>pscl</td><td>hurdle+summary</td><td>13.16</td><td>9.93x</td><td>1.37e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>1.43</td><td>MASS</td><td>glm.nb+summary</td><td>19.00</td><td>13.25x</td><td>4.01e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.08</td><td>stats</td><td>glm.fit+Wald</td><td>1.08</td><td>13.07x</td><td>3.4e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.09</td><td>stats</td><td>glm.fit+Wald(quasi)</td><td>1.32</td><td>14.68x</td><td>1.45e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.12</td><td>sandwich</td><td>glm+vcovHC</td><td>3.96</td><td>32.06x</td><td>4.69e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>3.18</td><td>pscl</td><td>zeroinfl(nb)+summary</td><td>40.20</td><td>12.62x</td><td>4.05e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>1.23</td><td>pscl</td><td>zeroinfl+summary</td><td>21.91</td><td>17.81x</td><td>1.25e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>1.81</td><td>betareg</td><td>betareg+summary</td><td>16.68</td><td>9.23x</td><td>1.81e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.10</td><td>stats</td><td>glm(quasi)+gcomp+Wald</td><td>2.70</td><td>26.68x</td><td>1.24e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.18</td><td>survival</td><td>coxph.fit(breslow)+Wald</td><td>0.47</td><td>2.61x</td><td>1.13e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>1.71</td><td>survival</td><td>coxph(null)+KM weighted residual mean diff + survdiff(rho=1)</td><td>5.25</td><td>3.06x</td><td>2.82e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>2.94</td><td>survival</td><td>survfit(median)+CI</td><td>3.22</td><td>1.1x</td><td>0.00984</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.01</td><td>survival</td><td>survdiff</td><td>1.90</td><td>139.62x</td><td>2.83e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.56</td><td>survival</td><td>coxph.fit(strat,breslow)+Wald</td><td>0.91</td><td>1.63x</td><td>0.0145</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.11</td><td>survival</td><td>survreg+summary</td><td>3.11</td><td>27.68x</td><td>1.05e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.17</td><td>VGAM</td><td>vglm+summary</td><td>14.87</td><td>87.4x</td><td>8.47e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.12</td><td>VGAM</td><td>vglm+summary</td><td>13.44</td><td>111.18x</td><td>7.27e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.69</td><td>ordinal</td><td>clm+gcomp+Wald</td><td>10.68</td><td>15.37x</td><td>0.000109</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>0.03</td><td>clinfun</td><td>jonckheere</td><td>0.47</td><td>14.91x</td><td>2.16e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.29</td><td>ordinal</td><td>clm+summary</td><td>6.24</td><td>21.16x</td><td>2.64e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>0.02</td><td>stats</td><td>mean(ridit)</td><td>0.28</td><td>15.04x</td><td>1e-05</td><td>***</td></tr>
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

EDI timings call exported C++ functions directly — no R6 objects are instantiated during benchmarking. As a result, **no R6 result caches exist to manage**. Each call to the C++ solver (e.g. `fast_logistic_regression_with_var_cpp`, `fast_ordinal_regression_cpp`) starts from a freshly zero-initialized parameter vector (or a model-specific data-driven initialization when `smart_cold_start = TRUE`). No prior-fit results are carried across timing repetitions, so every replication is a genuine cold start for the numerical optimizer.

<style>
    body, .markdown-body, .container {
        max-width: 1200px !important;
        width: 100% !important;
        margin: 0 auto !important;
    }
</style>
