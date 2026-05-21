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
    <tr style="background-color: #fff4bf;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.33</td><td>stats</td><td>t.test(pool)</td><td>0.46</td><td>1.37x</td><td>0.215</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>1.08</td><td>stats</td><td>wilcox.test</td><td>2.80</td><td>2.6x</td><td>0.000333</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.67</td><td>stats</td><td>lm.fit(interact)</td><td>1.34</td><td>2x</td><td>0.000262</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceContinOLS</td><td>continuous</td><td>0.40</td><td>stats</td><td>lm.fit</td><td>0.21</td><td>0.53x</td><td>0.0955</td><td></td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceContinQuantileRegr</td><td>continuous</td><td>2.50</td><td>quantreg</td><td>rq.fit</td><td>2.00</td><td>0.8x</td><td>0.191</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>1.25</td><td>MASS</td><td>rlm</td><td>3.29</td><td>2.63x</td><td>0.000283</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>1.83</td><td>stats</td><td>glm.fit(ident)</td><td>31.42</td><td>17.14x</td><td>6.8e-07</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceIncidExactFisher</td><td>incidence</td><td>2.00</td><td>stats</td><td>fisher.test</td><td>1.78</td><td>0.89x</td><td>0.365</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>1.22</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>3.39</td><td>2.78x</td><td>7.62e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>2.00</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>3.68</td><td>1.84x</td><td>0.00142</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>2.25</td><td>stats</td><td>glm.fit(log)</td><td>27.50</td><td>12.22x</td><td>2.17e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>2.00</td><td>stats</td><td>glm.fit</td><td>4.07</td><td>2.04x</td><td>6.42e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>1.00</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>1.68</td><td>1.68x</td><td>0.00322</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>1.67</td><td>stats</td><td>glm.fit(modified)</td><td>5.51</td><td>3.31x</td><td>2.72e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.10</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>1.61</td><td>16.13x</td><td>0.000213</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>6.38</td><td>stats</td><td>glm.fit(probit)</td><td>3.36</td><td>0.53x</td><td>0.0261</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.89</td><td>stats</td><td>prop.test</td><td>2.05</td><td>2.31x</td><td>0.00044</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>9.50</td><td>pscl</td><td>hurdle(nb)</td><td>90.67</td><td>9.54x</td><td>2.24e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>7.00</td><td>pscl</td><td>hurdle</td><td>50.00</td><td>7.14x</td><td>4.28e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>2.17</td><td>MASS</td><td>glm.nb</td><td>101.75</td><td>46.96x</td><td>2.29e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>3.00</td><td>stats</td><td>glm.fit</td><td>4.56</td><td>1.52x</td><td>0.00521</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>1.33</td><td>stats</td><td>glm.fit(quasi)</td><td>4.48</td><td>3.36x</td><td>0.000335</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>1.50</td><td>stats</td><td>glm.fit</td><td>4.56</td><td>3.04x</td><td>0.000299</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>10.15</td><td>pscl</td><td>zeroinfl(nb)</td><td>344.00</td><td>33.9x</td><td>1.02e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>7.68</td><td>pscl</td><td>zeroinfl</td><td>205.00</td><td>26.7x</td><td>1.47e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>24.12</td><td>betareg</td><td>betareg.fit</td><td>73.67</td><td>3.05x</td><td>0.000522</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>4.04</td><td>stats</td><td>glm.fit(quasi)</td><td>2.96</td><td>0.73x</td><td>0.04</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>1.50</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>3.18</td><td>2.12x</td><td>0.00109</td><td>**</td></tr>
    <tr style="background-color: #eceff1;"><td>InferencePropZeroOneInflatedBetaRegr</td><td>proportion</td><td>6.75</td><td>None</td><td>None</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>11.00</td><td>survival</td><td>coxph.fit</td><td>1.31</td><td>0.12x</td><td>5.78e-06</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceSurvivalDepCensTransformRegr</td><td>survival</td><td>22.75</td><td>None</td><td>None</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>10.00</td><td>survival</td><td>coxph(null)+KM weighted residual mean diff</td><td>12.68</td><td>1.27x</td><td>0.0147</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.33</td><td>survival</td><td>survfit(median)</td><td>8.32</td><td>24.95x</td><td>7.1e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.33</td><td>survival</td><td>survdiff</td><td>6.27</td><td>18.8x</td><td>1.57e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.12</td><td>survival</td><td>survfit(rmean)</td><td>6.18</td><td>49.61x</td><td>5.21e-07</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>112.50</td><td>survival</td><td>coxph.fit(strat)</td><td>1.17</td><td>0.01x</td><td>1.17e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>1.33</td><td>survival</td><td>survreg</td><td>8.85</td><td>6.64x</td><td>0.000339</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>2.75</td><td>VGAM</td><td>vglm(acat)</td><td>56.17</td><td>20.42x</td><td>2.54e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>10.75</td><td>ordinal</td><td>clm(cauchit)</td><td>42.64</td><td>3.97x</td><td>0.000233</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>11.50</td><td>ordinal</td><td>clm(cll)</td><td>24.87</td><td>2.16x</td><td>0.000386</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>2.50</td><td>VGAM</td><td>vglm(cratio)</td><td>49.90</td><td>19.96x</td><td>1.71e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>5.50</td><td>ordinal</td><td>clm+gcomp</td><td>32.83</td><td>5.97x</td><td>3.1e-05</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>74.50</td><td>clinfun</td><td>jonckheere</td><td>3.00</td><td>0.04x</td><td>3.78e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>11.00</td><td>ordinal</td><td>clm(probit)</td><td>18.13</td><td>1.65x</td><td>0.000386</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>7.00</td><td>ordinal</td><td>clm</td><td>23.33</td><td>3.33x</td><td>7.86e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>0.13</td><td>stats</td><td>mean(ridit)</td><td>0.63</td><td>5.06x</td><td>6.43e-07</td><td>***</td></tr>
  </tbody>
</table>


## Wald Test Performance (Full Inference)

This table compares the performance of **Full Inference** (Model Fit + Standard Error calculation + P-value derivation).
Unlike the point-estimation table above, these results include the computational cost of the variance-covariance matrix (Hessian or Fisher Information) and the Wald test statistic calculation.
All paths (EDI and Canonical) use a reduced sample size ($N=200$) for this full-inference benchmark to ensure iterative stability.
EDI timings in this table correspond to fixed `iBCRD` design objects.
EDI regression models (Logistic, Poisson) are benchmarked using the **IRLS** optimizer for these Wald tests.
**Note on Coverage**: `InferenceIncidCMH` is retained for coverage, but under `iBCRD` its EDI asymptotic p-value may be non-finite, in which case the row is reported as `NA`.
**Note on Slowdowns**: For some non-parametric tests (e.g. Jonckheere-Terpstra), EDI computes an **exact** p-value while R's counterpart uses a normal approximation for $N=200$, leading to a speedup < 1x.
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
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.18</td><td>stats</td><td>t.test(pool)</td><td>0.35</td><td>1.92x</td><td>0.000231</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.16</td><td>stats</td><td>wilcox.test</td><td>1.56</td><td>9.92x</td><td>1.3e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.62</td><td>stats</td><td>lm.fit(interact)+Wald</td><td>1.07</td><td>1.73x</td><td>0.00819</td><td>**</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceContinOLS</td><td>continuous</td><td>0.34</td><td>stats</td><td>lm.fit+Wald</td><td>0.16</td><td>0.46x</td><td>1.89e-05</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceContinQuantileRegr</td><td>continuous</td><td>5.31</td><td>quantreg</td><td>rq+summary</td><td>4.70</td><td>0.88x</td><td>0.0163</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.67</td><td>MASS</td><td>rlm+summary</td><td>3.13</td><td>4.67x</td><td>7.89e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidCMH</td><td>incidence</td><td>0.18</td><td>stats</td><td>mantelhaen</td><td>2.61</td><td>14.29x</td><td>2.17e-07</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidExactFisher</td><td>incidence</td><td>3.56</td><td>stats</td><td>fisher.test</td><td>1.71</td><td>0.48x</td><td>0.0013</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>1.84</td><td>stats</td><td>glm+gcomp(RD)+Wald</td><td>4.03</td><td>2.19x</td><td>3.37e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>1.75</td><td>stats</td><td>glm+gcomp(RR)+Wald</td><td>4.73</td><td>2.7x</td><td>9.97e-05</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>NA</td><td>stats</td><td>glm.fit+Wald(log)</td><td>6.38</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>3.09</td><td>stats</td><td>glm.fit+Wald</td><td>1.45</td><td>0.47x</td><td>0.000499</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.20</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>1.21</td><td>6.11x</td><td>4.04e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.46</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>1.53</td><td>3.32x</td><td>1.24e-05</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.81</td><td>stats</td><td>prop.test</td><td>1.12</td><td>1.38x</td><td>0.209</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>5.50</td><td>pscl</td><td>hurdle(nb)+summary</td><td>27.67</td><td>5.03x</td><td>3.24e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>4.88</td><td>pscl</td><td>hurdle+summary</td><td>18.50</td><td>3.79x</td><td>0.00129</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>2.68</td><td>MASS</td><td>glm.nb+summary</td><td>11.78</td><td>4.4x</td><td>5.79e-06</td><td>***</td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceCountPoisson</td><td>count</td><td>1.41</td><td>stats</td><td>glm.fit+Wald</td><td>1.83</td><td>1.3x</td><td>0.121</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.92</td><td>stats</td><td>glm.fit+Wald(quasi)</td><td>1.63</td><td>1.78x</td><td>0.00198</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>1.20</td><td>sandwich</td><td>glm+vcovHC</td><td>9.18</td><td>7.65x</td><td>0.000311</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>4.81</td><td>pscl</td><td>zeroinfl(nb)+summary</td><td>88.50</td><td>18.4x</td><td>9.79e-08</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>NA</td><td>pscl</td><td>zeroinfl+summary</td><td>45.56</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>4.95</td><td>betareg</td><td>betareg.fit+Wald</td><td>52.75</td><td>10.66x</td><td>1.32e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>2.22</td><td>stats</td><td>glm(quasi)+gcomp+Wald</td><td>4.37</td><td>1.97x</td><td>0.000352</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferencePropZeroOneInflatedBetaRegr</td><td>proportion</td><td>NA</td><td>None</td><td>None</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>106.75</td><td>survival</td><td>coxph.fit+Wald</td><td>1.04</td><td>0.01x</td><td>0.000121</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceSurvivalDepCensTransformRegr</td><td>survival</td><td>9.82</td><td>None</td><td>None</td><td>NA</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #fff4bf;"><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>13.67</td><td>survival</td><td>coxph(null)+KM weighted residual mean diff + survdiff(rho=1)</td><td>11.98</td><td>0.88x</td><td>0.246</td><td></td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>9.12</td><td>survival</td><td>survfit(median)+CI</td><td>6.47</td><td>0.71x</td><td>0.0421</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>1.65</td><td>survival</td><td>survdiff</td><td>3.87</td><td>2.34x</td><td>1.89e-05</td><td>***</td></tr>
    <tr style="background-color: #ffd9d9;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>29.50</td><td>survival</td><td>coxph.fit(strat)+Wald</td><td>1.01</td><td>0.03x</td><td>8.68e-06</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>1.66</td><td>survival</td><td>survreg+summary</td><td>6.84</td><td>4.13x</td><td>1.14e-07</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>1.92</td><td>VGAM</td><td>vglm+summary</td><td>30.00</td><td>15.66x</td><td>4.36e-06</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>NA</td><td>VGAM</td><td>vglm+summary</td><td>32.62</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>2.42</td><td>ordinal</td><td>clm+gcomp+Wald</td><td>14.29</td><td>5.92x</td><td>0.000132</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>NA</td><td>clinfun</td><td>jonckheere</td><td>0.87</td><td>NA</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>3.34</td><td>ordinal</td><td>clm+summary</td><td>10.60</td><td>3.17x</td><td>1.6e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>0.14</td><td>stats</td><td>mean(ridit)</td><td>0.37</td><td>2.63x</td><td>4.22e-07</td><td>***</td></tr>
  </tbody>
</table>

<style>
    body, .markdown-body, .container {
        max-width: 1200px !important;
        width: 100% !important;
        margin: 0 auto !important;
    }
</style>
