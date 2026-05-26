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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>1.34</td><td>stats</td><td>HL median pairwise diff</td><td>2.84</td><td>2.12x</td><td>1.38e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.07</td><td>stats</td><td>lm.fit</td><td>0.17</td><td>2.51x</td><td>1.3e-14</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceContinQuantileRegr</td><td>continuous</td><td>3.00</td><td>quantreg</td><td>rq.fit</td><td>2.52</td><td>0.84x</td><td>0.363</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.35</td><td>MASS</td><td>rlm(MM)</td><td>57.38</td><td>164.49x</td><td>2.53e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.53</td><td>stats</td><td>glm.fit(ident)</td><td>56.60</td><td>107.18x</td><td>1.15e-07</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>NA</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>1.82</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>NA</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>1.83</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>2.20</td><td>stats</td><td>glm.fit(log)</td><td>9.30</td><td>4.22x</td><td>1.52e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.20</td><td>stats</td><td>glm.fit</td><td>2.18</td><td>10.92x</td><td>1.44e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.24</td><td>stats</td><td>glm.fit(modified)</td><td>2.77</td><td>11.32x</td><td>4.92e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.84</td><td>stats</td><td>glm.fit(probit)</td><td>3.27</td><td>3.9x</td><td>1.07e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.07</td><td>stats</td><td>lm.fit(LPM)</td><td>0.18</td><td>2.75x</td><td>2.62e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>10.46</td><td>pscl</td><td>hurdle(nb)</td><td>59.12</td><td>5.65x</td><td>2.09e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>2.60</td><td>pscl</td><td>hurdle</td><td>29.00</td><td>11.14x</td><td>2.33e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>8.31</td><td>MASS</td><td>glm.nb</td><td>72.38</td><td>8.71x</td><td>7.65e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.21</td><td>stats</td><td>glm.fit</td><td>2.89</td><td>13.74x</td><td>4.08e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.19</td><td>stats</td><td>glm.fit(quasi)</td><td>2.46</td><td>12.66x</td><td>5.12e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.19</td><td>stats</td><td>glm.fit</td><td>2.53</td><td>13.54x</td><td>4.24e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>2.17</td><td>pscl</td><td>zeroinfl(nb)</td><td>206.00</td><td>94.81x</td><td>1.17e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>5.95</td><td>pscl</td><td>zeroinfl</td><td>101.83</td><td>17.13x</td><td>3.27e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>7.56</td><td>betareg</td><td>betareg.fit</td><td>45.25</td><td>5.98x</td><td>4.41e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.19</td><td>stats</td><td>glm.fit(quasi)</td><td>2.28</td><td>12.07x</td><td>3.1e-08</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>NA</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>2.45</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.26</td><td>survival</td><td>coxph.fit(breslow)</td><td>0.35</td><td>1.34x</td><td>5.70e-39</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.21</td><td>survival</td><td>survfit(median)</td><td>5.21</td><td>24.77x</td><td>3.1e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.07</td><td>survival</td><td>survdiff</td><td>2.91</td><td>43.88x</td><td>5.29e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.21</td><td>survival</td><td>survfit(rmean)</td><td>3.51</td><td>16.43x</td><td>5.48e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.33</td><td>survival</td><td>coxph.fit(strat)</td><td>0.45</td><td>1.35x</td><td>9.74e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.11</td><td>survival</td><td>survreg</td><td>5.54</td><td>48.73x</td><td>5.15e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>2.05</td><td>VGAM</td><td>vglm(acat)</td><td>138.00</td><td>67.39x</td><td>5.99e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>1.09</td><td>ordinal</td><td>clm(cauchit)</td><td>14.28</td><td>13.09x</td><td>1.22e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>1.38</td><td>ordinal</td><td>clm(cloglog)</td><td>11.95</td><td>8.65x</td><td>1.26e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>2.36</td><td>VGAM</td><td>vglm(cratio)</td><td>29.00</td><td>12.27x</td><td>2.62e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>1.10</td><td>ordinal</td><td>clm+gcomp</td><td>19.29</td><td>17.53x</td><td>4.53e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>1.28</td><td>ordinal</td><td>clm(probit)</td><td>11.65</td><td>9.1x</td><td>9.47e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>1.08</td><td>ordinal</td><td>clm</td><td>11.09</td><td>10.31x</td><td>3.65e-11</td><td>***</td></tr>
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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.10</td><td>stats</td><td>t.test(pool)</td><td>0.16</td><td>1.54x</td><td>1.74e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.14</td><td>stats</td><td>wilcox.test</td><td>0.64</td><td>4.67x</td><td>1.88e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.32</td><td>stats</td><td>lm.fit(interact)+Wald</td><td>0.59</td><td>1.85x</td><td>5.84e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.03</td><td>stats</td><td>lm.fit+Wald</td><td>0.13</td><td>3.65x</td><td>3.75e-91</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceContinQuantileRegr</td><td>continuous</td><td>3.20</td><td>quantreg</td><td>rq+summary</td><td>2.77</td><td>0.87x</td><td>0.061</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.41</td><td>MASS</td><td>rlm+summary</td><td>1.69</td><td>4.09x</td><td>4.28e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidCMH</td><td>incidence</td><td>0.10</td><td>stats</td><td>mantelhaen</td><td>1.24</td><td>11.93x</td><td>1.19e-09</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidExactFisher</td><td>incidence</td><td>1.43</td><td>stats</td><td>fisher.test</td><td>0.82</td><td>0.57x</td><td>2.08e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.70</td><td>stats</td><td>glm+gcomp(RD)+Wald</td><td>2.38</td><td>3.38x</td><td>3.8e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.70</td><td>stats</td><td>glm+gcomp(RR)+Wald</td><td>2.27</td><td>3.23x</td><td>1.5e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>0.85</td><td>stats</td><td>glm.fit+Wald(log)</td><td>1.79</td><td>2.1x</td><td>1.02e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.53</td><td>stats</td><td>glm.fit+Wald</td><td>0.69</td><td>1.31x</td><td>1.43e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.11</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>0.67</td><td>5.9x</td><td>4.39e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.17</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>0.84</td><td>5.01x</td><td>1.51e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.24</td><td>stats</td><td>prop.test</td><td>0.41</td><td>1.69x</td><td>5.61e-13</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>24.00</td><td>pscl</td><td>hurdle(nb)+summary</td><td>24.13</td><td>1.01x</td><td>0.573</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>3.24</td><td>pscl</td><td>hurdle+summary</td><td>13.46</td><td>4.16x</td><td>6.26e-06</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceCountNegBin</td><td>count</td><td>567.50</td><td>MASS</td><td>glm.nb+summary</td><td>6.05</td><td>0.01x</td><td>2.9e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.84</td><td>stats</td><td>glm.fit+Wald</td><td>0.88</td><td>1.06x</td><td>0.0234</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.70</td><td>stats</td><td>glm.fit+Wald(quasi)</td><td>0.95</td><td>1.35x</td><td>1.68e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.55</td><td>sandwich</td><td>glm+vcovHC</td><td>4.15</td><td>7.5x</td><td>2.03e-14</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>NA</td><td>pscl</td><td>zeroinfl(nb)+summary</td><td>80.17</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>2.70</td><td>pscl</td><td>zeroinfl+summary</td><td>54.88</td><td>20.34x</td><td>5.82e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>2.12</td><td>betareg</td><td>betareg.fit+Wald</td><td>29.17</td><td>13.78x</td><td>9.8e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.85</td><td>stats</td><td>glm(quasi)+gcomp+Wald</td><td>2.19</td><td>2.56x</td><td>5.09e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.31</td><td>survival</td><td>coxph.fit(breslow)+Wald</td><td>0.39</td><td>1.27x</td><td>2.06e-21</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>5.53</td><td>survival</td><td>coxph(null)+KM weighted residual mean diff + survdiff(rho=1)</td><td>5.05</td><td>0.91x</td><td>0.00927</td><td>**</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>4.63</td><td>survival</td><td>survfit(median)+CI</td><td>3.43</td><td>0.74x</td><td>1e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>1.14</td><td>survival</td><td>survdiff</td><td>1.81</td><td>1.59x</td><td>4.9e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.36</td><td>survival</td><td>coxph.fit(strat,breslow)+Wald</td><td>0.42</td><td>1.19x</td><td>1.60e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.69</td><td>survival</td><td>survreg+summary</td><td>3.76</td><td>5.43x</td><td>8.9e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.90</td><td>VGAM</td><td>vglm+summary</td><td>18.94</td><td>21.04x</td><td>2.96e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.76</td><td>VGAM</td><td>vglm+summary</td><td>15.46</td><td>20.32x</td><td>9.16e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.91</td><td>ordinal</td><td>clm+gcomp+Wald</td><td>7.62</td><td>8.4x</td><td>1.48e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>0.19</td><td>clinfun</td><td>jonckheere</td><td>0.47</td><td>2.5x</td><td>1.92e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.86</td><td>ordinal</td><td>clm+summary</td><td>6.88</td><td>8.04x</td><td>1.11e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>0.14</td><td>stats</td><td>mean(ridit)</td><td>0.16</td><td>1.15x</td><td>0.0016</td><td>**</td></tr>
  </tbody>
</table>

<style>
    body, .markdown-body, .container {
        max-width: 1200px !important;
        width: 100% !important;
        margin: 0 auto !important;
    }
</style>
