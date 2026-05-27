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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>1.19</td><td>stats</td><td>HL median pairwise diff</td><td>2.26</td><td>1.9x</td><td>1.57e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.06</td><td>stats</td><td>lm.fit</td><td>0.15</td><td>2.31x</td><td>1.34e-08</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceContinQuantileRegr</td><td>continuous</td><td>3.00</td><td>quantreg</td><td>rq.fit</td><td>2.64</td><td>0.88x</td><td>0.438</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.36</td><td>MASS</td><td>rlm(MM)</td><td>59.13</td><td>166.13x</td><td>6.76e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.24</td><td>stats</td><td>glm.fit(ident)</td><td>17.68</td><td>75.06x</td><td>5.77e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.28</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>1.77</td><td>6.43x</td><td>2.51e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.30</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>1.83</td><td>6.12x</td><td>7.74e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>2.23</td><td>stats</td><td>glm.fit(log)</td><td>8.05</td><td>3.61x</td><td>1.73e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.18</td><td>stats</td><td>glm.fit</td><td>1.91</td><td>10.73x</td><td>3.55e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.21</td><td>stats</td><td>glm.fit(modified)</td><td>2.55</td><td>12.23x</td><td>3.58e-11</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>NA</td><td>stats</td><td>glm.fit(probit)</td><td>2.86</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.07</td><td>stats</td><td>lm.fit(LPM)</td><td>0.18</td><td>2.64x</td><td>8.58e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>11.96</td><td>pscl</td><td>hurdle(nb)</td><td>66.17</td><td>5.53x</td><td>3.69e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>2.69</td><td>pscl</td><td>hurdle</td><td>25.71</td><td>9.56x</td><td>9.43e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>8.08</td><td>MASS</td><td>glm.nb</td><td>61.62</td><td>7.63x</td><td>1.25e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.20</td><td>stats</td><td>glm.fit</td><td>2.35</td><td>12.04x</td><td>3.23e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.20</td><td>stats</td><td>glm.fit(quasi)</td><td>2.41</td><td>12.03x</td><td>1.73e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.20</td><td>stats</td><td>glm.fit</td><td>2.75</td><td>14.07x</td><td>4.95e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>2.79</td><td>pscl</td><td>zeroinfl(nb)</td><td>209.50</td><td>75.17x</td><td>4.53e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>5.95</td><td>pscl</td><td>zeroinfl</td><td>98.25</td><td>16.51x</td><td>2.26e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>7.19</td><td>betareg</td><td>betareg.fit</td><td>35.92</td><td>4.99x</td><td>1.05e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.21</td><td>stats</td><td>glm.fit(quasi)</td><td>2.03</td><td>9.53x</td><td>2.06e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.31</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>2.43</td><td>7.87x</td><td>8.72e-08</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>21.25</td><td>survival</td><td>coxph.fit(breslow)</td><td>0.72</td><td>0.03x</td><td>1.71e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.21</td><td>survival</td><td>survfit(median)</td><td>5.02</td><td>24x</td><td>4.02e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.06</td><td>survival</td><td>survdiff</td><td>2.67</td><td>42.86x</td><td>1.43e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.22</td><td>survival</td><td>survfit(rmean)</td><td>3.59</td><td>16.69x</td><td>1.35e-09</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>17.75</td><td>survival</td><td>coxph.fit(strat)</td><td>0.52</td><td>0.03x</td><td>8.84e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.11</td><td>survival</td><td>survreg</td><td>5.47</td><td>47.74x</td><td>1.49e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.81</td><td>VGAM</td><td>vglm(acat)</td><td>41.12</td><td>50.67x</td><td>8.67e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>0.95</td><td>ordinal</td><td>clm(cauchit)</td><td>15.87</td><td>16.64x</td><td>2.52e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>1.38</td><td>ordinal</td><td>clm(cloglog)</td><td>10.87</td><td>7.88x</td><td>1.34e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.80</td><td>VGAM</td><td>vglm(cratio)</td><td>28.57</td><td>35.85x</td><td>4.18e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>1.28</td><td>ordinal</td><td>clm+gcomp</td><td>18.64</td><td>14.58x</td><td>8.85e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>1.37</td><td>ordinal</td><td>clm(probit)</td><td>10.86</td><td>7.9x</td><td>2.41e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>1.07</td><td>ordinal</td><td>clm</td><td>9.68</td><td>9.06x</td><td>1.86e-12</td><td>***</td></tr>
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
**Row Highlighting**: Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light red rows indicate `Speedup < 1` and `Timing Pval < 0.05`; light yellow rows indicate `Timing Pval >= 0.05`; light grey rows indicate `NA` timing comparisons.

<table>
  <thead>
    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.16</td><td>stats</td><td>t.test(pool)</td><td>0.23</td><td>1.43x</td><td>1.9e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.26</td><td>stats</td><td>wilcox.test</td><td>0.95</td><td>3.67x</td><td>6.74e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>1.52</td><td>stats</td><td>lm.fit(interact)+Wald</td><td>2.69</td><td>1.78x</td><td>0.00807</td><td>**</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceContinOLS</td><td>continuous</td><td>0.19</td><td>stats</td><td>lm.fit+Wald</td><td>0.08</td><td>0.43x</td><td>4.73e-09</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceContinQuantileRegr</td><td>continuous</td><td>3.74</td><td>quantreg</td><td>rq+summary</td><td>2.98</td><td>0.8x</td><td>4.37e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.45</td><td>MASS</td><td>rlm+summary</td><td>1.90</td><td>4.19x</td><td>3.6e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidCMH</td><td>incidence</td><td>0.16</td><td>stats</td><td>mantelhaen</td><td>1.95</td><td>12.56x</td><td>5.46e-10</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidExactFisher</td><td>incidence</td><td>1.38</td><td>stats</td><td>fisher.test</td><td>1.27</td><td>0.92x</td><td>0.013</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>1.81</td><td>stats</td><td>glm+gcomp(RD)+Wald</td><td>10.52</td><td>5.8x</td><td>4.09e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>3.85</td><td>stats</td><td>glm+gcomp(RR)+Wald</td><td>7.61</td><td>1.98x</td><td>0.00242</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>0.73</td><td>stats</td><td>glm.fit+Wald(log)</td><td>2.11</td><td>2.88x</td><td>4.94e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.60</td><td>stats</td><td>glm.fit+Wald</td><td>0.73</td><td>1.21x</td><td>0.000482</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.46</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>2.83</td><td>6.19x</td><td>1.01e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.86</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>2.82</td><td>3.28x</td><td>2.87e-05</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>1.29</td><td>stats</td><td>prop.test</td><td>1.41</td><td>1.1x</td><td>0.0856</td><td></td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>66.33</td><td>pscl</td><td>hurdle(nb)+summary</td><td>24.21</td><td>0.36x</td><td>2.43e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>3.71</td><td>pscl</td><td>hurdle+summary</td><td>13.13</td><td>3.54x</td><td>2.03e-06</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceCountNegBin</td><td>count</td><td>630.00</td><td>MASS</td><td>glm.nb+summary</td><td>6.71</td><td>0.01x</td><td>8.52e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.85</td><td>stats</td><td>glm.fit+Wald</td><td>0.98</td><td>1.15x</td><td>0.000291</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.91</td><td>stats</td><td>glm.fit+Wald(quasi)</td><td>1.64</td><td>1.81x</td><td>0.00507</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>3.18</td><td>sandwich</td><td>glm+vcovHC</td><td>21.65</td><td>6.8x</td><td>6.73e-05</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>NA</td><td>pscl</td><td>zeroinfl(nb)+summary</td><td>43.58</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>2.86</td><td>pscl</td><td>zeroinfl+summary</td><td>65.67</td><td>22.93x</td><td>4.51e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>2.47</td><td>betareg</td><td>betareg.fit+Wald</td><td>27.72</td><td>11.24x</td><td>8.16e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>3.76</td><td>stats</td><td>glm(quasi)+gcomp+Wald</td><td>9.75</td><td>2.59x</td><td>0.000112</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>9.28</td><td>survival</td><td>coxph.fit(breslow)+Wald</td><td>0.53</td><td>0.06x</td><td>2.23e-11</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>7.50</td><td>survival</td><td>coxph(null)+KM weighted residual mean diff + survdiff(rho=1)</td><td>8.21</td><td>1.09x</td><td>0.465</td><td></td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>20.13</td><td>survival</td><td>survfit(median)+CI</td><td>11.85</td><td>0.59x</td><td>0.0109</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>1.39</td><td>survival</td><td>survdiff</td><td>2.84</td><td>2.04x</td><td>7.69e-09</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>13.89</td><td>survival</td><td>coxph.fit(strat,breslow)+Wald</td><td>1.29</td><td>0.09x</td><td>1.02e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>2.06</td><td>survival</td><td>survreg+summary</td><td>4.78</td><td>2.32x</td><td>3.44e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.85</td><td>VGAM</td><td>vglm+summary</td><td>20.06</td><td>23.71x</td><td>7.65e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.83</td><td>VGAM</td><td>vglm+summary</td><td>18.46</td><td>22.24x</td><td>1.77e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>4.14</td><td>ordinal</td><td>clm+gcomp+Wald</td><td>29.42</td><td>7.1x</td><td>5.68e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>0.25</td><td>clinfun</td><td>jonckheere</td><td>0.60</td><td>2.42x</td><td>1.35e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.94</td><td>ordinal</td><td>clm+summary</td><td>5.85</td><td>6.21x</td><td>1.43e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>0.53</td><td>stats</td><td>mean(ridit)</td><td>0.76</td><td>1.43x</td><td>0.0147</td><td>*</td></tr>
  </tbody>
</table>

<style>
    body, .markdown-body, .container {
        max-width: 1200px !important;
        width: 100% !important;
        margin: 0 auto !important;
    }
</style>
