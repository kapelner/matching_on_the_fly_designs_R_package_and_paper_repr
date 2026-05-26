# EDI Exhaustive C++ Model Fit Benchmarks

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

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
*   **Row Highlighting:** Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light red rows indicate `Speedup < 1` and `Timing Pval < 0.05`; light yellow rows indicate `Timing Pval >= 0.05`; light grey rows indicate `NA` timing comparisons.
*   **Constraints**: Matched-pair/KK and highly custom paths are excluded as per user request.

## Results

<table>
  <thead>
    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr style="background-color: #fff4bf;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>3.36</td><td>stats</td><td>HL median pairwise diff</td><td>3.14</td><td>0.93x</td><td>0.126</td><td></td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceContinOLS</td><td>continuous</td><td>0.43</td><td>stats</td><td>lm.fit</td><td>0.17</td><td>0.4x</td><td>1.57e-10</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceContinQuantileRegr</td><td>continuous</td><td>5.00</td><td>quantreg</td><td>rq.fit</td><td>2.82</td><td>0.56x</td><td>0.00228</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>10.79</td><td>MASS</td><td>rlm(MM)</td><td>64.17</td><td>5.94x</td><td>1.12e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>6.15</td><td>stats</td><td>glm.fit(ident)</td><td>35.40</td><td>5.76x</td><td>1.72e-09</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>5.49</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>2.95</td><td>0.54x</td><td>0.00229</td><td>**</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>6.15</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>3.09</td><td>0.5x</td><td>6.79e-11</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>173.50</td><td>stats</td><td>glm.fit(log)</td><td>17.34</td><td>0.1x</td><td>4.5e-09</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>4.50</td><td>stats</td><td>glm.fit</td><td>2.39</td><td>0.53x</td><td>5.98e-12</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>7.12</td><td>stats</td><td>glm.fit(modified)</td><td>3.09</td><td>0.43x</td><td>2e-11</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>10.68</td><td>stats</td><td>glm.fit(probit)</td><td>3.38</td><td>0.32x</td><td>3.32e-08</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.47</td><td>stats</td><td>lm.fit(LPM)</td><td>0.20</td><td>0.43x</td><td>1.88e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>34.93</td><td>pscl</td><td>hurdle(nb)</td><td>81.00</td><td>2.32x</td><td>6.79e-07</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>51.25</td><td>pscl</td><td>hurdle</td><td>30.92</td><td>0.6x</td><td>2.07e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>26.94</td><td>MASS</td><td>glm.nb</td><td>72.83</td><td>2.7x</td><td>2.65e-06</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceCountPoisson</td><td>count</td><td>8.70</td><td>stats</td><td>glm.fit</td><td>3.21</td><td>0.37x</td><td>1.1e-13</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>11.29</td><td>stats</td><td>glm.fit(quasi)</td><td>3.12</td><td>0.28x</td><td>4.8e-10</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceCountRobustPoisson</td><td>count</td><td>9.97</td><td>stats</td><td>glm.fit</td><td>3.60</td><td>0.36x</td><td>2.98e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>71.50</td><td>pscl</td><td>zeroinfl(nb)</td><td>331.50</td><td>4.64x</td><td>3.72e-05</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>198.75</td><td>pscl</td><td>zeroinfl</td><td>120.83</td><td>0.61x</td><td>0.000127</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>19.55</td><td>betareg</td><td>betareg.fit</td><td>41.00</td><td>2.1x</td><td>1.48e-11</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>6.03</td><td>stats</td><td>glm.fit(quasi)</td><td>3.68</td><td>0.61x</td><td>1.88e-06</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>5.33</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>2.64</td><td>0.5x</td><td>2.43e-09</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>2.67</td><td>survival</td><td>coxph.fit(breslow)</td><td>1.01</td><td>0.38x</td><td>0.000177</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.24</td><td>survival</td><td>survfit(median)</td><td>7.00</td><td>29.17x</td><td>4.19e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.30</td><td>survival</td><td>survdiff</td><td>3.79</td><td>12.44x</td><td>1.39e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.22</td><td>survival</td><td>survfit(rmean)</td><td>4.49</td><td>20.17x</td><td>1.78e-08</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>2.00</td><td>survival</td><td>coxph.fit(strat)</td><td>0.57</td><td>0.29x</td><td>4.38e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>1.85</td><td>survival</td><td>survreg</td><td>6.80</td><td>3.68x</td><td>2.09e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>46.60</td><td>VGAM</td><td>vglm(acat)</td><td>71.83</td><td>1.54x</td><td>1.22e-07</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>34.00</td><td>ordinal</td><td>clm(cauchit)</td><td>19.18</td><td>0.56x</td><td>8.16e-11</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>26.29</td><td>ordinal</td><td>clm(cloglog)</td><td>18.86</td><td>0.72x</td><td>4.44e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>12.86</td><td>VGAM</td><td>vglm(cratio)</td><td>35.20</td><td>2.74x</td><td>4.25e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>30.86</td><td>ordinal</td><td>clm+gcomp</td><td>38.88</td><td>1.26x</td><td>0.011</td><td>*</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>28.93</td><td>ordinal</td><td>clm(probit)</td><td>13.08</td><td>0.45x</td><td>3.37e-13</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>22.27</td><td>ordinal</td><td>clm</td><td>10.82</td><td>0.49x</td><td>3.72e-12</td><td>***</td></tr>
  </tbody>
</table>

<style>
    body, .markdown-body, .container {
        max-width: 1200px !important;
        width: 100% !important;
        margin: 0 auto !important;
    }
</style>
