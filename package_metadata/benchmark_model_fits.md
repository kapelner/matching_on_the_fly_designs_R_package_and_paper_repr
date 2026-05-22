# EDI Exhaustive C++ Model Fit Benchmarks

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Benchmark Dataset Specification

All benchmarks were performed on a synthetic clinical-trial-scale dataset generated for each response type. The data generation process ensures numerical stability and fair solver comparison by using the following parameters:

*   **Sample Size ($N$):** 1,000 subjects for most models; 500 subjects for survival models. Exact and trend tests may use smaller scaled samples (N=100-500) as noted in the results.
*   **Predictors ($p$):** 5 total predictors, including a global intercept, a balanced binary treatment assignment from fixed `iBCRD`, and 4 continuous covariates ($X \sim \text{Normal}(0, 1)$).
*   **Effect Sizes:** Covariate coefficients are sampled from $\text{Normal}(0, 0.5)$, matching the warm-start benchmark data model. The treatment coefficient is set to 0.5 in the linear predictor so the benchmarked treatment effect is meaningfully separated from zero.
*   **EDI Design Template:** EDI benchmark objects are instantiated on a fixed `iBCRD` design.
*   **Response Generation:**
    *   **Continuous:** Linear model with additive $\text{Normal}(0, 0.5)$ noise.
    *   **Incidence:** Binary outcomes via a Logistic link.
    *   **Count:** Integer outcomes via Poisson or Negative Binomial distributions with an exponential link.
    *   **Proportion:** Continuous outcomes in $(0, 1)$ via a Beta distribution with a logit link.
    *   **Survival:** Exponentially distributed event times with approximately 20% random censoring.
    *   **Ordinal:** 3-level categorical outcomes generated from the same ordinal construction used in the warm-start benchmark.

## Methodology

*   **Pure Solver Timing:** Results reflect the time taken for the core numerical optimization. We exclude R6 object instantiation, design matrix construction, and standard error estimation (which often uses different R-side matrix inversion logic) to isolate the efficiency of the underlying C++ backends.
*   **Solver-Only Prebuilds:** For EDI rows, benchmark setup prebuilds exposed observed-data design matrices, reduced design matrices, and other fixed working inputs outside the timed region when the implementation exposes those hooks. Canonical rows using low-level matrix interfaces are likewise timed on prebuilt inputs.
*   **Smart Cold Starts:** EDI models were initialized with `smart_cold_start = TRUE`, utilizing package-optimized heuristic starting values.
*   **Randomization Design:** EDI timings in this table correspond to `iBCRD` design objects.
*   **Low-Level Comparison:** Canonical R timings use the fastest available internal interfaces (e.g., `.fit` functions) to remove R's formula parsing and environment management overhead.
*   **Limitation:** Some canonical comparators only expose formula-based APIs rather than comparable low-level fit kernels. Those rows remain included, but their canonical timings may still contain formula/model-frame overhead beyond the numerical solver itself.
*   **Averaging:** All timings are medians over 10 warmed runs measured with adaptive batched `system.time`; paths below 0.01 ms use `microbenchmark(times = 2000)` instead.
*   **Timing P-Value:** `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.
*   **Row Highlighting:** Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light red rows indicate `Speedup < 1` and `Timing Pval < 0.05`; light yellow rows indicate `Timing Pval >= 0.05`; light grey rows indicate `NA` timing comparisons.
*   **Constraints**: Matched-pair/KK and highly custom paths are excluded as per user request.

## Results

<table>
  <thead>
    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.05</td><td>stats</td><td>t.test(pool)</td><td>0.17</td><td>3.25x</td><td>5.11e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.46</td><td>stats</td><td>HL median pairwise diff</td><td>3.47</td><td>7.57x</td><td>2.29e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.29</td><td>stats</td><td>lm.fit(interact)</td><td>0.65</td><td>2.26x</td><td>6.53e-11</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceContinOLS</td><td>continuous</td><td>0.30</td><td>stats</td><td>lm.fit</td><td>0.15</td><td>0.51x</td><td>0.00372</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinQuantileRegr</td><td>continuous</td><td>1.40</td><td>quantreg</td><td>rq.fit</td><td>1.49</td><td>1.06x</td><td>0.0291</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.50</td><td>MASS</td><td>rlm(MM)</td><td>42.90</td><td>85.8x</td><td>7.04e-16</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.69</td><td>stats</td><td>glm.fit(ident)</td><td>13.74</td><td>19.98x</td><td>7.63e-19</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceIncidExactFisher</td><td>incidence</td><td>0.64</td><td>stats</td><td>fisher.test</td><td>0.68</td><td>1.06x</td><td>0.532</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.60</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>1.56</td><td>2.59x</td><td>7.53e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.67</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>1.61</td><td>2.41x</td><td>0.000338</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>0.88</td><td>stats</td><td>glm.fit(log)</td><td>16.63</td><td>19.01x</td><td>3.93e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>1.00</td><td>stats</td><td>glm.fit</td><td>2.00</td><td>2x</td><td>1.3e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.17</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>0.66</td><td>3.93x</td><td>1.32e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.62</td><td>stats</td><td>glm.fit(modified)</td><td>2.23</td><td>3.57x</td><td>7.64e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.04</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>0.79</td><td>19.89x</td><td>9.29e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.60</td><td>stats</td><td>glm.fit(probit)</td><td>2.24</td><td>3.73x</td><td>1.12e-21</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.41</td><td>stats</td><td>lm.fit(LPM)</td><td>0.15</td><td>0.37x</td><td>7.23e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>2.17</td><td>pscl</td><td>hurdle(nb)</td><td>47.20</td><td>21.78x</td><td>4.76e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>2.67</td><td>pscl</td><td>hurdle</td><td>21.88</td><td>8.2x</td><td>1.05e-23</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>1.17</td><td>MASS</td><td>glm.nb</td><td>56.88</td><td>48.75x</td><td>8.55e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>1.40</td><td>stats</td><td>glm.fit</td><td>2.31</td><td>1.65x</td><td>7.83e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.58</td><td>stats</td><td>glm.fit(quasi)</td><td>2.07</td><td>3.54x</td><td>1.46e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.67</td><td>stats</td><td>glm.fit</td><td>2.07</td><td>3.1x</td><td>0.000183</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>3.40</td><td>pscl</td><td>zeroinfl(nb)</td><td>164.75</td><td>48.52x</td><td>1.47e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>2.32</td><td>pscl</td><td>zeroinfl</td><td>77.17</td><td>33.3x</td><td>1.31e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>10.29</td><td>betareg</td><td>betareg.fit</td><td>32.00</td><td>3.11x</td><td>4.57e-18</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.75</td><td>stats</td><td>glm.fit(quasi)</td><td>2.05</td><td>2.74x</td><td>7.23e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.67</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>1.64</td><td>2.46x</td><td>6.73e-09</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.60</td><td>survival</td><td>coxph.fit(breslow)</td><td>0.69</td><td>1.15x</td><td>0.363</td><td></td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>5.33</td><td>survival</td><td>coxph(null)+KM weighted residual mean diff</td><td>4.38</td><td>0.82x</td><td>0.00114</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.14</td><td>survival</td><td>survfit(median)</td><td>3.74</td><td>26.19x</td><td>9.47e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.11</td><td>survival</td><td>survdiff</td><td>2.09</td><td>18.8x</td><td>7.59e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.05</td><td>survival</td><td>survfit(rmean)</td><td>3.26</td><td>61.43x</td><td>1.65e-08</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.67</td><td>survival</td><td>coxph.fit(strat)</td><td>0.45</td><td>0.68x</td><td>0.000132</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.72</td><td>survival</td><td>survreg</td><td>4.35</td><td>6.02x</td><td>7.98e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>1.91</td><td>VGAM</td><td>vglm(acat)</td><td>31.17</td><td>16.33x</td><td>1.1e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>3.93</td><td>ordinal</td><td>clm(cauchit)</td><td>11.06</td><td>2.81x</td><td>0.0002</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>3.43</td><td>ordinal</td><td>clm(cloglog)</td><td>9.00</td><td>2.62x</td><td>9.95e-18</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>1.22</td><td>VGAM</td><td>vglm(cratio)</td><td>21.00</td><td>17.18x</td><td>2.27e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.67</td><td>ordinal</td><td>clm+gcomp</td><td>15.42</td><td>23.13x</td><td>4.32e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>3.21</td><td>ordinal</td><td>clm(probit)</td><td>8.23</td><td>2.56x</td><td>8.14e-18</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>2.57</td><td>ordinal</td><td>clm</td><td>9.13</td><td>3.55x</td><td>4.17e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>0.09</td><td>stats</td><td>mean(ridit; ref=control)</td><td>0.23</td><td>2.56x</td><td>4.33e-07</td><td>***</td></tr>
  </tbody>
</table>


## Wald Test Performance (Full Inference)

This table compares the performance of **Full Inference** (Model Fit + Standard Error calculation + P-value derivation).
Unlike the point-estimation table above, these results include the computational cost of the variance-covariance matrix (Hessian or Fisher Information) and the Wald test statistic calculation.
All paths (EDI and Canonical) use a reduced sample size ($N=200$) for this full-inference benchmark to ensure iterative stability.
EDI timings in this table correspond to fixed `iBCRD` design objects.
EDI regression models (Logistic, Poisson) are benchmarked using the **IRLS** optimizer for these Wald tests.
**Note on Coverage**: `InferenceIncidCMH` is retained for coverage, but under `iBCRD` its EDI asymptotic p-value may be non-finite, in which case the row is reported as `NA`.
**Note on Accessors**: EDI asymptotic Wald timings in this table explicitly call both `compute_asymp_confidence_interval()` and `compute_asymp_two_sided_pval()` so the benchmark includes the full asymptotic inference accessor path.
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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.11</td><td>stats</td><td>t.test(pool)</td><td>0.14</td><td>1.24x</td><td>1.05e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.08</td><td>stats</td><td>wilcox.test</td><td>0.59</td><td>7.05x</td><td>2e-17</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.31</td><td>stats</td><td>lm.fit(interact)+Wald</td><td>0.53</td><td>1.69x</td><td>6.13e-07</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceContinOLS</td><td>continuous</td><td>0.21</td><td>stats</td><td>lm.fit+Wald</td><td>0.07</td><td>0.35x</td><td>1.29e-18</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceContinQuantileRegr</td><td>continuous</td><td>2.52</td><td>quantreg</td><td>rq+summary</td><td>1.89</td><td>0.75x</td><td>0.00114</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.30</td><td>MASS</td><td>rlm+summary</td><td>1.47</td><td>4.93x</td><td>9.88e-16</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidCMH</td><td>incidence</td><td>0.12</td><td>stats</td><td>mantelhaen</td><td>1.14</td><td>9.87x</td><td>3.34e-15</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidExactFisher</td><td>incidence</td><td>1.36</td><td>stats</td><td>fisher.test</td><td>0.76</td><td>0.56x</td><td>6.27e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.60</td><td>stats</td><td>glm+gcomp(RD)+Wald</td><td>1.98</td><td>3.3x</td><td>7.47e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.62</td><td>stats</td><td>glm+gcomp(RR)+Wald</td><td>2.01</td><td>3.26x</td><td>1.96e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>0.61</td><td>stats</td><td>glm.fit+Wald(log)</td><td>1.59</td><td>2.58x</td><td>3.78e-16</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.93</td><td>stats</td><td>glm.fit+Wald</td><td>0.60</td><td>0.64x</td><td>1.79e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.16</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>0.63</td><td>3.94x</td><td>1.13e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.20</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>0.75</td><td>3.69x</td><td>3.73e-12</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.41</td><td>stats</td><td>prop.test</td><td>0.41</td><td>1x</td><td>0.458</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>7.06</td><td>pscl</td><td>hurdle(nb)+summary</td><td>15.36</td><td>2.18x</td><td>4.57e-15</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>1.99</td><td>pscl</td><td>hurdle+summary</td><td>8.35</td><td>4.19x</td><td>4.21e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>1.16</td><td>MASS</td><td>glm.nb+summary</td><td>5.72</td><td>4.95x</td><td>1.54e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.61</td><td>stats</td><td>glm.fit+Wald</td><td>0.78</td><td>1.28x</td><td>1.8e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.43</td><td>stats</td><td>glm.fit+Wald(quasi)</td><td>0.81</td><td>1.87x</td><td>1.38e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.51</td><td>sandwich</td><td>glm+vcovHC</td><td>3.74</td><td>7.38x</td><td>6.19e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>2.38</td><td>pscl</td><td>zeroinfl(nb)+summary</td><td>68.83</td><td>28.89x</td><td>7.27e-13</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>2.26</td><td>pscl</td><td>zeroinfl+summary</td><td>20.05</td><td>8.88x</td><td>1.55e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>2.29</td><td>betareg</td><td>betareg.fit+Wald</td><td>21.85</td><td>9.55x</td><td>6.96e-16</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.79</td><td>stats</td><td>glm(quasi)+gcomp+Wald</td><td>2.09</td><td>2.65x</td><td>8.23e-13</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.87</td><td>survival</td><td>coxph.fit(breslow)+Wald</td><td>0.43</td><td>0.49x</td><td>2.94e-15</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>4.76</td><td>survival</td><td>coxph(null)+KM weighted residual mean diff + survdiff(rho=1)</td><td>4.53</td><td>0.95x</td><td>0.000493</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>4.02</td><td>survival</td><td>survfit(median)+CI</td><td>3.01</td><td>0.75x</td><td>7.67e-10</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.70</td><td>survival</td><td>survdiff</td><td>1.50</td><td>2.13x</td><td>2.91e-14</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.49</td><td>survival</td><td>coxph.fit(strat,breslow)+Wald</td><td>0.34</td><td>0.7x</td><td>1.92e-09</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.70</td><td>survival</td><td>survreg+summary</td><td>2.74</td><td>3.9x</td><td>9.88e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.82</td><td>VGAM</td><td>vglm+summary</td><td>12.00</td><td>14.71x</td><td>3.1e-16</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.74</td><td>VGAM</td><td>vglm+summary</td><td>12.62</td><td>17.02x</td><td>2.54e-14</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.85</td><td>ordinal</td><td>clm+gcomp+Wald</td><td>7.00</td><td>8.19x</td><td>8.5e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>0.15</td><td>clinfun</td><td>jonckheere</td><td>0.39</td><td>2.54x</td><td>7.63e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>1.45</td><td>ordinal</td><td>clm+summary</td><td>4.70</td><td>3.25x</td><td>1.09e-11</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>0.09</td><td>stats</td><td>mean(ridit)</td><td>0.17</td><td>1.87x</td><td>2.88e-16</td><td>***</td></tr>
  </tbody>
</table>

<style>
    body, .markdown-body, .container {
        max-width: 1200px !important;
        width: 100% !important;
        margin: 0 auto !important;
    }
</style>
