# EDI Exhaustive C++ Model Fit Benchmarks

_Generated: 2026-06-10 09:10:35 EDT_

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Compilation Context

These rows are read from build metadata compiled into the loaded `EDI` shared object via `edi_build_info_cpp()`.

**Compilation warning:** EDI model-fit timings are sensitive to the compiler flags used to build the loaded `EDI.so`. If EDI is compiled without the proper optimized flags, or with flags that are known to degrade these kernels such as problematic LTO builds, the benchmark can show substantial performance regressions that reflect the binary build rather than the modeling algorithms.

*   **EDI shared object:** `/home/kapelner/R/x86_64-pc-linux-gnu-library/4.7/EDI/libs/EDI.so`
*   **EDI shared object mtime:** `2026-06-09 01:21:43`
*   **Capture method:** `configure-generated header compiled into EDI.so`
*   **Build timestamp:** `2026-06-05 09:55:23 EDT`
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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.84</td><td>stats</td><td>HL median pairwise diff</td><td>5.73</td><td>6.86x</td><td>0.000159</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.05</td><td>stats</td><td>lm.fit</td><td>0.27</td><td>5.67x</td><td>1.51e-07</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>1.73</td><td>quantreg</td><td>rq.fit</td><td>1.87</td><td>1.08x</td><td>0.583</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.42</td><td>MASS</td><td>rlm(MM)</td><td>67.38</td><td>161.14x</td><td>1.64e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.20</td><td>stats</td><td>glm.fit(ident)</td><td>17.94</td><td>91.59x</td><td>1.82e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.38</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>4.79</td><td>12.6x</td><td>6.32e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.37</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>3.72</td><td>9.95x</td><td>2.91e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>5.20</td><td>stats</td><td>glm.fit(log)</td><td>13.03</td><td>2.51x</td><td>0.000456</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.26</td><td>stats</td><td>glm.fit</td><td>4.20</td><td>15.93x</td><td>6.11e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.38</td><td>stats</td><td>glm.fit(modified)</td><td>4.91</td><td>12.9x</td><td>1.48e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.86</td><td>stats</td><td>glm.fit(probit)</td><td>3.00</td><td>3.5x</td><td>5.15e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.05</td><td>stats</td><td>lm.fit(LPM)</td><td>0.27</td><td>5.72x</td><td>1.01e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>14.46</td><td>pscl</td><td>hurdle(nb)</td><td>86.00</td><td>5.95x</td><td>5.93e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>3.56</td><td>pscl</td><td>hurdle</td><td>42.90</td><td>12.07x</td><td>0.00012</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>10.87</td><td>MASS</td><td>glm.nb</td><td>100.75</td><td>9.27x</td><td>3.51e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.42</td><td>stats</td><td>glm.fit</td><td>4.16</td><td>9.92x</td><td>1.32e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.43</td><td>stats</td><td>glm.fit(quasi)</td><td>3.37</td><td>7.87x</td><td>2.42e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.44</td><td>stats</td><td>glm.fit</td><td>4.18</td><td>9.47x</td><td>1.17e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>4.06</td><td>pscl</td><td>zeroinfl(nb)</td><td>370.50</td><td>91.36x</td><td>1.34e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>9.60</td><td>pscl</td><td>zeroinfl</td><td>169.50</td><td>17.66x</td><td>1.39e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>10.26</td><td>betareg</td><td>betareg.fit</td><td>61.00</td><td>5.94x</td><td>1.19e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.23</td><td>stats</td><td>glm.fit(quasi)</td><td>2.66</td><td>11.65x</td><td>1.65e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.35</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>3.58</td><td>10.31x</td><td>0.00606</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.68</td><td>survival</td><td>coxph.fit(breslow)</td><td>1.58</td><td>2.32x</td><td>2.59e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.03</td><td>survival</td><td>survfit(median)</td><td>8.87</td><td>316.8x</td><td>8.86e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.04</td><td>survival</td><td>survdiff</td><td>5.17</td><td>138.52x</td><td>1.97e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.03</td><td>survival</td><td>survfit(rmean)</td><td>7.68</td><td>235.52x</td><td>1.19e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>1.59</td><td>survival</td><td>coxph.fit(strat)</td><td>1.96</td><td>1.23x</td><td>0.00858</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.16</td><td>survival</td><td>survreg</td><td>7.42</td><td>45.72x</td><td>6.47e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.73</td><td>VGAM</td><td>vglm(acat)</td><td>18.75</td><td>25.74x</td><td>1.97e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>0.94</td><td>ordinal</td><td>clm(cauchit)</td><td>14.27</td><td>15.24x</td><td>3.2e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>0.94</td><td>ordinal</td><td>clm(cloglog)</td><td>11.94</td><td>12.74x</td><td>2.05e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.61</td><td>VGAM</td><td>vglm(cratio)</td><td>18.64</td><td>30.48x</td><td>4.52e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>1.51</td><td>ordinal</td><td>clm+gcomp</td><td>35.92</td><td>23.71x</td><td>7.74e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>1.03</td><td>ordinal</td><td>clm(probit)</td><td>11.00</td><td>10.7x</td><td>7.84e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>1.46</td><td>ordinal</td><td>clm</td><td>18.44</td><td>12.61x</td><td>9.25e-08</td><td>***</td></tr>
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
