# EDI Exhaustive C++ Model Fit Benchmarks

_Generated: 2026-06-22 10:14:38 HST_

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Compilation Context

These rows are read from build metadata compiled into the loaded `EDI` shared object via `edi_build_info_cpp()`.

**Compilation warning:** EDI model-fit timings are sensitive to the compiler flags used to build the loaded `EDI.so`. If EDI is compiled without the proper optimized flags, or with flags that are known to degrade these kernels such as problematic LTO builds, the benchmark can show substantial performance regressions that reflect the binary build rather than the modeling algorithms.

*   **EDI shared object:** `/home/kapelner/R/x86_64-pc-linux-gnu-library/4.7/EDI/libs/EDI.so`
*   **EDI shared object mtime:** `2026-06-22 09:43:01`
*   **Capture method:** `configure-generated header compiled into EDI.so`
*   **Build timestamp:** `2026-06-22 08:15:15 HST`
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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.64</td><td>stats</td><td>HL median pairwise diff</td><td>2.01</td><td>3.12x</td><td>7.67e-23</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.03</td><td>stats</td><td>lm.fit</td><td>0.15</td><td>5.47x</td><td>1.39e-28</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>1.69</td><td>quantreg</td><td>rq.fit</td><td>1.68</td><td>1x</td><td>0.554</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.30</td><td>MASS</td><td>rlm(MM)</td><td>43.00</td><td>143.07x</td><td>2.92e-19</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.15</td><td>stats</td><td>glm.fit(ident)</td><td>13.38</td><td>90.47x</td><td>5.11e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.20</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>1.76</td><td>8.95x</td><td>3.46e-23</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.18</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>1.74</td><td>9.79x</td><td>5.72e-23</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>2.11</td><td>stats</td><td>glm.fit(log)</td><td>4.97</td><td>2.36x</td><td>1.35e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.14</td><td>stats</td><td>glm.fit</td><td>1.94</td><td>13.99x</td><td>5.37e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.26</td><td>stats</td><td>glm.fit(modified)</td><td>2.02</td><td>7.71x</td><td>1.52e-23</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.57</td><td>stats</td><td>glm.fit(probit)</td><td>2.02</td><td>3.56x</td><td>1.79e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.03</td><td>stats</td><td>lm.fit(LPM)</td><td>0.14</td><td>5.44x</td><td>3.63e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>9.35</td><td>pscl</td><td>hurdle(nb)</td><td>51.20</td><td>5.48x</td><td>2.43e-31</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>2.16</td><td>pscl</td><td>hurdle</td><td>23.32</td><td>10.81x</td><td>9.96e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>6.25</td><td>MASS</td><td>glm.nb</td><td>55.50</td><td>8.88x</td><td>1.37e-31</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.24</td><td>stats</td><td>glm.fit</td><td>2.43</td><td>10.08x</td><td>7.09e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.25</td><td>stats</td><td>glm.fit(quasi)</td><td>2.33</td><td>9.32x</td><td>2.82e-23</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.29</td><td>stats</td><td>glm.fit</td><td>1.79</td><td>6.23x</td><td>5.26e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>2.34</td><td>pscl</td><td>zeroinfl(nb)</td><td>182.75</td><td>78.1x</td><td>4.75e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>5.35</td><td>pscl</td><td>zeroinfl</td><td>83.83</td><td>15.67x</td><td>3.08e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>6.62</td><td>betareg</td><td>betareg.fit</td><td>36.50</td><td>5.52x</td><td>2.21e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.15</td><td>stats</td><td>glm.fit(quasi)</td><td>1.50</td><td>10.12x</td><td>1.4e-22</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.18</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>1.57</td><td>8.98x</td><td>2.19e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.40</td><td>survival</td><td>coxph.fit(breslow)</td><td>0.69</td><td>1.72x</td><td>8.25e-29</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.02</td><td>survival</td><td>survfit(median)</td><td>3.98</td><td>263.4x</td><td>2.93e-30</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.02</td><td>survival</td><td>survdiff</td><td>2.04</td><td>88.76x</td><td>2.51e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.02</td><td>survival</td><td>survfit(rmean)</td><td>2.66</td><td>175.36x</td><td>1.07e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.69</td><td>survival</td><td>coxph.fit(strat)</td><td>0.76</td><td>1.09x</td><td>3.69e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.11</td><td>survival</td><td>survreg</td><td>4.44</td><td>40.58x</td><td>8.3e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.59</td><td>VGAM</td><td>vglm(acat)</td><td>14.90</td><td>25.12x</td><td>1.22e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>0.80</td><td>ordinal</td><td>clm(cauchit)</td><td>10.14</td><td>12.68x</td><td>3.46e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>0.73</td><td>ordinal</td><td>clm(cloglog)</td><td>8.87</td><td>12.21x</td><td>9.24e-30</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.49</td><td>VGAM</td><td>vglm(cratio)</td><td>14.00</td><td>28.69x</td><td>1.69e-35</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.78</td><td>ordinal</td><td>clm+gcomp</td><td>14.30</td><td>18.25x</td><td>8.08e-36</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>0.87</td><td>ordinal</td><td>clm(probit)</td><td>7.80</td><td>8.91x</td><td>2.59e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.73</td><td>ordinal</td><td>clm</td><td>9.00</td><td>12.28x</td><td>7.04e-33</td><td>***</td></tr>
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
