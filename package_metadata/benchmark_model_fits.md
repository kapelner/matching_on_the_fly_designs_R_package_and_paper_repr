# EDI Exhaustive C++ Model Fit Benchmarks

_Generated: 2026-05-29 01:34:02 EDT_

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Compilation Context

These rows are read from build metadata compiled into the loaded `EDI` shared object via `edi_build_info_cpp()`.

**Compilation warning:** EDI model-fit timings are sensitive to the compiler flags used to build the loaded `EDI.so`. If EDI is compiled without the proper optimized flags, or with flags that are known to degrade these kernels such as problematic LTO builds, the benchmark can show substantial performance regressions that reflect the binary build rather than the modeling algorithms.

*   **EDI shared object:** `/home/kapelner/R/x86_64-pc-linux-gnu-library/4.7/EDI/libs/EDI.so`
*   **EDI shared object mtime:** `2026-05-29 00:26:35`
*   **Capture method:** `configure-generated header compiled into EDI.so`
*   **Build timestamp:** `2026-05-29 00:17:19 EDT`
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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.56</td><td>stats</td><td>HL median pairwise diff</td><td>2.59</td><td>4.6x</td><td>4.71e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.01</td><td>stats</td><td>lm.fit</td><td>0.14</td><td>9.52x</td><td>2.15e-11</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>1.50</td><td>quantreg</td><td>rq.fit</td><td>1.31</td><td>0.87x</td><td>0.991</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.17</td><td>MASS</td><td>rlm(MM)</td><td>40.60</td><td>237.34x</td><td>1.08e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.19</td><td>stats</td><td>glm.fit(ident)</td><td>12.55</td><td>67.8x</td><td>4.99e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.18</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>1.83</td><td>10.32x</td><td>1.3e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.16</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>1.93</td><td>12.17x</td><td>1.01e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>0.13</td><td>stats</td><td>glm.fit(log)</td><td>15.12</td><td>118.47x</td><td>1.37e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.13</td><td>stats</td><td>glm.fit</td><td>1.67</td><td>12.95x</td><td>1.7e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.15</td><td>stats</td><td>glm.fit(modified)</td><td>2.18</td><td>14.99x</td><td>5.48e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.63</td><td>stats</td><td>glm.fit(probit)</td><td>2.08</td><td>3.31x</td><td>1.44e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.02</td><td>stats</td><td>lm.fit(LPM)</td><td>0.18</td><td>11.22x</td><td>2.09e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>8.80</td><td>pscl</td><td>hurdle(nb)</td><td>45.50</td><td>5.17x</td><td>2.68e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>1.96</td><td>pscl</td><td>hurdle</td><td>20.50</td><td>10.47x</td><td>1.15e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>6.18</td><td>MASS</td><td>glm.nb</td><td>49.75</td><td>8.05x</td><td>1.39e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.15</td><td>stats</td><td>glm.fit</td><td>1.99</td><td>12.91x</td><td>3.27e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.14</td><td>stats</td><td>glm.fit(quasi)</td><td>1.79</td><td>12.38x</td><td>2.89e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.18</td><td>stats</td><td>glm.fit</td><td>2.10</td><td>11.84x</td><td>1.36e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>2.05</td><td>pscl</td><td>zeroinfl(nb)</td><td>159.50</td><td>77.8x</td><td>1.8e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>4.83</td><td>pscl</td><td>zeroinfl</td><td>72.83</td><td>15.07x</td><td>5.25e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>5.66</td><td>betareg</td><td>betareg.fit</td><td>33.40</td><td>5.9x</td><td>2.1e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.17</td><td>stats</td><td>glm.fit(quasi)</td><td>1.92</td><td>11.57x</td><td>4.33e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.16</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>1.83</td><td>11.53x</td><td>4.04e-07</td><td>***</td></tr>
    <tr><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.67</td><td>survival</td><td>coxph.fit(breslow)</td><td>0.67</td><td>1.01x</td><td>0.748</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.02</td><td>survival</td><td>survfit(median)</td><td>4.07</td><td>196.77x</td><td>3.97e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.02</td><td>survival</td><td>survdiff</td><td>2.51</td><td>102.2x</td><td>3.61e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.02</td><td>survival</td><td>survfit(rmean)</td><td>3.04</td><td>168.31x</td><td>5e-10</td><td>***</td></tr>
    <tr><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.50</td><td>survival</td><td>coxph.fit(strat)</td><td>0.44</td><td>0.88x</td><td>0.506</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.06</td><td>survival</td><td>survreg</td><td>3.92</td><td>70.23x</td><td>7.11e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.58</td><td>VGAM</td><td>vglm(acat)</td><td>29.33</td><td>50.67x</td><td>1.47e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>0.77</td><td>ordinal</td><td>clm(cauchit)</td><td>11.79</td><td>15.26x</td><td>4.46e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>1.11</td><td>ordinal</td><td>clm(cloglog)</td><td>9.07</td><td>8.14x</td><td>6.48e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.40</td><td>VGAM</td><td>vglm(cratio)</td><td>20.91</td><td>52.42x</td><td>8.99e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.86</td><td>ordinal</td><td>clm+gcomp</td><td>16.50</td><td>19.2x</td><td>1.69e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>0.96</td><td>ordinal</td><td>clm(probit)</td><td>8.54</td><td>8.89x</td><td>4.88e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.71</td><td>ordinal</td><td>clm</td><td>7.65</td><td>10.71x</td><td>1.57e-11</td><td>***</td></tr>
  </tbody>
</table>

## Garbage Collection and Cache Management

To ensure that the benchmark results are highly precise, reproducible, and represent the actual computation speed of the numerical solvers, the benchmarking harness uses the following garbage collection and cache management strategies:

### 1. Garbage Collection (GC) Filtering
Garbage collection cycles run automatically by the R interpreter and can introduce significant, arbitrary pauses that skew timing measurements. To isolate the execution time of the code from R's GC overhead:
* **GC Disabling**: We disable R's memory stress-testing mode using `gctorture(FALSE)` before running timing loops.
* **Proactive Compaction**: In the `system.time()` path, we invoke `gc(verbose = FALSE)` immediately before timing each replicate. This starts the timer on a clean, compacted heap, minimizing the likelihood of triggering an automatic garbage collection cycle mid-replicate.
* **Automatic Filtering**: In the microbenchmarking path, we utilize the `bench::mark()` engine with the `filter_gc = TRUE` parameter, which automatically tracks and discards timing iterations during which a garbage collection event occurred.

### 2. Fair Cache Management for R6 Estimators
Many EDI estimators utilize R6 objects that cache model fits (`cached_mod`) and computed estimates/variances (`cached_values`) to avoid redundant computations on subsequent accessor calls.
* **Clean Calculations**: If these caches were not cleared between benchmark repetitions, iterations 2 to $N$ would return the cached results instantly in $O(1)$ time, which would prevent measuring the actual numerical optimization speed.
* **Cache Cleansing**: To compare the raw C++ optimization algorithms against R's canonical solvers fairly (e.g. in simulation or bootstrap loops where the data/weights change on every run and the cache is not reusable), the R6 estimator caches are explicitly cleared on every single repetition.
* **Overhead Subtraction**: Modifying environments and setting variables to `NULL` in R introduces small computational overhead. To ensure this cleanup cost is not counted against EDI, the benchmarking harness isolates the cleanup time by timing it separately (via `setup_only` control runs) and subtracting it from the total execution time, ensuring a clean measurement of the solver itself.


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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.10</td><td>stats</td><td>t.test(pool)</td><td>0.16</td><td>1.58x</td><td>6.27e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.11</td><td>stats</td><td>wilcox.test</td><td>0.97</td><td>9.12x</td><td>1.12e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.35</td><td>stats</td><td>lm.fit(interact)+Wald</td><td>0.48</td><td>1.38x</td><td>4.5e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.05</td><td>stats</td><td>lm.fit+Wald</td><td>0.08</td><td>1.52x</td><td>2.97e-05</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>2.65</td><td>quantreg</td><td>rq+summary</td><td>2.09</td><td>0.79x</td><td>0.0105</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.30</td><td>MASS</td><td>rlm+summary</td><td>1.36</td><td>4.57x</td><td>2.43e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidCMH</td><td>incidence</td><td>0.10</td><td>stats</td><td>mantelhaen</td><td>1.12</td><td>11.72x</td><td>1.54e-10</td><td>***</td></tr>
    <tr><td>InferenceIncidExactFisher</td><td>incidence</td><td>0.99</td><td>stats</td><td>fisher.test</td><td>0.81</td><td>0.82x</td><td>0.000984</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.55</td><td>stats</td><td>glm+gcomp(RD)+Wald</td><td>1.85</td><td>3.36x</td><td>6.07e-17</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.49</td><td>stats</td><td>glm+gcomp(RR)+Wald</td><td>1.89</td><td>3.85x</td><td>8.32e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>1.18</td><td>stats</td><td>glm.fit+Wald(log)</td><td>6.11</td><td>5.18x</td><td>7.9e-17</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.26</td><td>stats</td><td>glm.fit+Wald</td><td>0.72</td><td>2.72x</td><td>1.47e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.10</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>0.68</td><td>6.53x</td><td>2.04e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.14</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>0.75</td><td>5.43x</td><td>1.49e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.53</td><td>stats</td><td>glm.fit(probit)+Wald</td><td>0.90</td><td>1.68x</td><td>3.61e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.20</td><td>stats</td><td>prop.test</td><td>0.41</td><td>2.07x</td><td>7.22e-08</td><td>***</td></tr>
    <tr><td>InferenceCountHurdleNegBin</td><td>count</td><td>18.54</td><td>pscl</td><td>hurdle(nb)+summary</td><td>13.83</td><td>0.75x</td><td>4.66e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>2.77</td><td>pscl</td><td>hurdle+summary</td><td>10.03</td><td>3.62x</td><td>4.36e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>3.18</td><td>MASS</td><td>glm.nb+summary</td><td>16.44</td><td>5.17x</td><td>3.51e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.55</td><td>stats</td><td>glm.fit+Wald</td><td>0.83</td><td>1.53x</td><td>4.28e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.57</td><td>stats</td><td>glm.fit+Wald(quasi)</td><td>0.93</td><td>1.62x</td><td>6.7e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.50</td><td>sandwich</td><td>glm+vcovHC</td><td>3.58</td><td>7.23x</td><td>3.04e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>5.34</td><td>pscl</td><td>zeroinfl(nb)+summary</td><td>30.37</td><td>5.69x</td><td>3.15e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>2.61</td><td>pscl</td><td>zeroinfl+summary</td><td>16.56</td><td>6.34x</td><td>6.79e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>1.98</td><td>betareg</td><td>betareg.fit+Wald</td><td>25.21</td><td>12.77x</td><td>5.02e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.61</td><td>stats</td><td>glm(quasi)+gcomp+Wald</td><td>1.90</td><td>3.12x</td><td>2.13e-10</td><td>***</td></tr>
    <tr><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.67</td><td>survival</td><td>coxph.fit(breslow)+Wald</td><td>0.44</td><td>0.66x</td><td>0.000302</td><td>***</td></tr>
    <tr><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>5.22</td><td>survival</td><td>coxph(null)+KM weighted residual mean diff + survdiff(rho=1)</td><td>5.12</td><td>0.98x</td><td>0.923</td><td></td></tr>
    <tr><td>InferenceSurvivalKMDiff</td><td>survival</td><td>2.98</td><td>survival</td><td>survfit(median)+CI</td><td>3.13</td><td>1.05x</td><td>0.0791</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.94</td><td>survival</td><td>survdiff</td><td>1.59</td><td>1.69x</td><td>2.83e-10</td><td>***</td></tr>
    <tr><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.35</td><td>survival</td><td>coxph.fit(strat,breslow)+Wald</td><td>0.33</td><td>0.94x</td><td>0.0321</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.46</td><td>survival</td><td>survreg+summary</td><td>3.06</td><td>6.64x</td><td>2.84e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.54</td><td>VGAM</td><td>vglm+summary</td><td>14.97</td><td>27.48x</td><td>1.35e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.52</td><td>VGAM</td><td>vglm+summary</td><td>13.85</td><td>26.44x</td><td>7.49e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.62</td><td>ordinal</td><td>clm+gcomp+Wald</td><td>6.21</td><td>10.07x</td><td>1.96e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>0.17</td><td>clinfun</td><td>jonckheere</td><td>0.45</td><td>2.56x</td><td>5.85e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.89</td><td>ordinal</td><td>clm+summary</td><td>5.47</td><td>6.15x</td><td>4.48e-17</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>0.06</td><td>stats</td><td>mean(ridit)</td><td>0.17</td><td>2.97x</td><td>1.07e-10</td><td>***</td></tr>
  </tbody>
</table>

## Garbage Collection and Cache Management

To ensure that the benchmark results are highly precise, reproducible, and represent the actual computation speed of the numerical solvers, the benchmarking harness uses the following garbage collection and cache management strategies:

### 1. Garbage Collection (GC) Filtering
Garbage collection cycles run automatically by the R interpreter and can introduce significant, arbitrary pauses that skew timing measurements. To isolate the execution time of the code from R's GC overhead:
* **GC Disabling**: We disable R's memory stress-testing mode using `gctorture(FALSE)` before running timing loops.
* **Proactive Compaction**: In the `system.time()` path, we invoke `gc(verbose = FALSE)` immediately before timing each replicate. This starts the timer on a clean, compacted heap, minimizing the likelihood of triggering an automatic garbage collection cycle mid-replicate.
* **Automatic Filtering**: In the microbenchmarking path, we utilize the `bench::mark()` engine with the `filter_gc = TRUE` parameter, which automatically tracks and discards timing iterations during which a garbage collection event occurred.

### 2. Fair Cache Management for R6 Estimators
Many EDI estimators utilize R6 objects that cache model fits (`cached_mod`) and computed estimates/variances (`cached_values`) to avoid redundant computations on subsequent accessor calls.
* **Clean Calculations**: If these caches were not cleared between benchmark repetitions, iterations 2 to $N$ would return the cached results instantly in $O(1)$ time, which would prevent measuring the actual numerical optimization speed.
* **Cache Cleansing**: To compare the raw C++ optimization algorithms against R's canonical solvers fairly (e.g. in simulation or bootstrap loops where the data/weights change on every run and the cache is not reusable), the R6 estimator caches are explicitly cleared on every single repetition.
* **Overhead Subtraction**: Modifying environments and setting variables to `NULL` in R introduces small computational overhead. To ensure this cleanup cost is not counted against EDI, the benchmarking harness isolates the cleanup time by timing it separately (via `setup_only` control runs) and subtracting it from the total execution time, ensuring a clean measurement of the solver itself.

<style>
    body, .markdown-body, .container {
        max-width: 1200px !important;
        width: 100% !important;
        margin: 0 auto !important;
    }
</style>
