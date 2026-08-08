# Smart Cold Start Strategies in EDI

This report summarizes how the "smart cold start" (initial parameter selection from no prior knowledge) is performed for each likelihood path in the `EDI` package. 

A **Cold Start** occurs when fitting a model on a new dataset with no previous iteration data. **Smart Starts** use optimized heuristics (like OLS) to find a starting point near the MLE, significantly reducing iterations compared to a naive start at zero.

## Heuristic Strategies by Model Path

Strategies marked with "(R Standard)" have been aligned with core R functions like `glm.fit` and `survival::survreg`.

| Model / Likelihood Path | Primary Strategy for $\beta$ (Treatment & Covariates) | Strategy for Secondary Parameters ($\alpha, \sigma, \theta$, etc.) |
| :--- | :--- | :--- |
| **Logistic Regression** | **OLS on $\text{logit}(\frac{y + 0.5}{2})$**: (R Standard) Shrinks 0/1 response toward 0.5 to provide a robust initial separation plane. | N/A |
| **Poisson Regression** | **OLS on $\ln(y + 0.1)$**: (R Standard) Uses a small offset for zero counts to allow log-linear initialization. | N/A |
| **Negative Binomial** | **OLS on $\ln(y + 1)$**: Log-linear approximation using OLS. | $\ln(\theta) = 0$ (initial dispersion set to 1). |
| **Weibull AFT** | **OLS on $\ln(y) + 0.572$**: (R
 Standard) Uses uncensored observations with a Gumbel mean shift. | **Moment-based scale**: $\ln(\sigma) = 0.5 \ln(\text{resid\_var} / 1.64)$. |
| **Beta Regression** | **OLS on $\text{logit}(y)$**: Maps (0,1) response to real line via logit. | $\ln(\phi) = 2.0$. |
| **Ordinal Regression** | **OLS on $y$**: Uses ordinal levels as a continuous response for $\beta$. | **Quantile Mapping**: $\alpha$ intercepts set by mapping sample quantiles. |

---

## 2026 Updated Benchmark: End-to-End Smart Start (N=500, Strong Signal)

This benchmark evaluates the performance of **Naive (Zero)** vs. **Smart** initialization for one-off estimations on high-signal data.

**Benchmark Environment:**
* **Sample Size ($n$):** 500
* **Covariates ($p$):** 5
* **Signal Strength:** Very Strong ($\beta \sim N(0, 2.0^2)$).
* **Heuristic Engine:** Internal C++ optimized OLS.

### Performance Summary (End-to-End)

| Model Family | Naive Time (ms) | Smart Time (ms) | Smart Its | Naive Its | Speedup (%) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Weibull (AFT)** | 0.32 ms | 0.16 ms | **5** | 24 | **+50.0%** |
| **Beta Regression** | 13.68 ms | 7.76 ms | **N/A** | N/A | **+43.3%** |
| **ZINB** | 15.88 ms | 10.20 ms | **N/A** | N/A | **+35.8%** |
| **Logistic (IRLS)** | 0.24 ms | 0.24 ms | 6 | 6 | **0.0%** |
| **Poisson (IRLS)** | 0.08 ms | 0.20 ms | 6 | **2** | *-150%* |

---

## Key Findings: Robustness vs. Speed

### 1. High-Complexity Acceleration
For models with complex likelihood surfaces or many nuisance parameters (e.g., **Weibull, Beta, ZINB**), a smart start is transformative. Providing a statistically motivated initial guess reduces the search space for the L-BFGS or Newton solver, typically cutting execution time by **35-50%**.

### 2. Boundary Safety
For binomial models (Log-Binomial, Identity Link), smart starts are an essential **insurance policy**. A naive start at zero often puts the initial gradient in a region where the mean parameter violates the $[0,1]$ boundary. The smart start ensures the optimizer begins in the feasible region.

### 3. Simple GLM Efficiency
For simple models with a **Strong Signal**, a naive start at zero can sometimes be superior. In the Poisson test above, the zero-start converged in 2 iterations, while the OLS-start landed further away and took 6. This confirms that for very fast models (<0.5ms), the overhead of smart-start logic can occasionally exceed its benefit.

**Conclusion:** EDI implements **Smart Starts** as the global default because the reliability gain (ensuring convergence on difficult data) far outweighs the sub-millisecond overhead in simple cases.

---

## Implementation TODOs

### Fix table rows that contradict the actual default policy

`EDI/R/globals.R:638-651` (`get_cold_start_dispatch_policy()`) overrides `smart_cold_start` to `FALSE` (naive zero start) by default for `InferenceIncidLogRegr` (Logistic Regression) and `InferenceCountPoisson` (Poisson Regression), specifically *because* benchmarking showed the OLS warm-up is net-negative for IRLS at typical trial sizes — this is exactly the "Poisson (IRLS)" `-150%` finding already in this document's benchmark table (line 41), just never connected back to the row above it.

- [ ] TODO-1: Rewrite the **Logistic Regression** and **Poisson Regression** table rows (lines 13-14) to state that the OLS heuristic described is available but **not the default** — the default cold start is a plain zero vector per `EDI/R/globals.R:638-651`. The heuristic itself is still real and implemented at `EDI/src/fast_logistic_regression.cpp:117` (`edi_opt::logistic_smart_cold_start(X, y)`) and `EDI/src/fast_poisson_regression.cpp:107,121`; it's just gated off by policy.
- [ ] TODO-2: Add a footnote/subsection listing every other class the same policy overrides to `FALSE` by default, so the table isn't misread as "smart start always on": `InferencePropFractionalLogit`, `InferencePropGComp*`, `InferenceIncidGComp*`, `InferenceIncidKKGComp*`, `InferenceCountQuasiPoisson`, `InferenceCountRobustPoisson`, `InferenceIncidModifiedPoisson` (full list at `EDI/R/globals.R:643-650`).
- [ ] TODO-3: Reconcile section "3. Simple GLM Efficiency" (line 53-54) with TODO-1 — it currently reads as an observation, not as the explanation for an already-shipped policy decision. Rewrite it to say the package resolved this by disabling smart-start for exactly these families, rather than leaving it as an open tension.

### Fix the missing ZINB row

- [ ] TODO-4: Add a **ZINB** row to the "Heuristic Strategies by Model Path" table (it's benchmarked at line 39 but has no row above). Per `EDI/src/fast_zinb.cpp:175-190`, ZINB's actual smart-start is **not** OLS-based like NegBin's: the count-model intercept is set to `log(mean(y))` (clamped at `1e-8`), all other count-model and zero-inflation coefficients start at `0`, and `log(theta) = 0`. Do not describe it as "OLS on ln(y+1)" — that's NegBin's heuristic, not ZINB's.

### Expand table coverage — heuristics already implemented but undocumented

Every one of these families already has a working `smart_cold_start` branch in its C++ kernel (default `TRUE` unless noted in TODO-2), but has no row in this document's table. Each needs its own row with the precise heuristic read out of the cited source, the same way the existing 6 rows do — do not guess the formula from the family name.

- [ ] TODO-5: **Hurdle NegBin** / **Hurdle Poisson (GLMM)** — `EDI/src/fast_hurdle_negbin.cpp` (search `smart_cold_start` branches near its `fast_hurdle_negbin_internal`), `EDI/src/fast_hurdle_poisson_glmm.cpp:490`, `EDI/src/fast_zero_augmented_poisson.cpp:286`.
- [ ] TODO-6: **Zero-One-Inflated Beta** (proportion family) — `EDI/src/fast_zero_one_inflated_beta.cpp:435-457` (uses `ols_smart_cold_start_beta` on the zero-mass and one-mass logit sub-models).
- [ ] TODO-7: **Log-Binomial** (incidence, `InferenceIncidLogBinom`) — `EDI/src/fast_log_binomial_regression.cpp:185-190,358-363` (two variants: `ols_smart_cold_start_beta_on_log1p_or_legacy` and `ols_smart_cold_start_beta_or_legacy` — determine which is the current default and why two exist).
- [ ] TODO-8: **Probit Regression** (incidence) — `EDI/src/fast_probit_regression.cpp:138`.
- [ ] TODO-9: **Robust Regression** (continuous, M-estimator) — `EDI/src/fast_robust_regression.cpp:101`.
- [ ] TODO-10: **Cox PH** and **Stratified Cox PH** (survival) — `EDI/src/fast_coxph_regression.cpp:267,433`.
- [ ] TODO-11: **General parametric survival models** (used by `InferenceSurvivalKKClaytonCopulaOneLik` and related KK survival paths) — `EDI/src/fast_survival_models_optim.cpp:681`.
- [ ] TODO-12: **Logistic GLMM** and **Poisson GLMM** (KK-design mixed models) — `EDI/src/fast_logistic_glmm.cpp:481`, `EDI/src/fast_poisson_glmm.cpp:380`.
- [ ] TODO-13: **Ordinal GLMM**, **Ordinal Cauchit**, **Ordinal Cloglog**, **Ordinal Ordered Probit**, **Adjacent-Category Logit**, **Continuation-Ratio Regression**, **Stereotype Logit** — these are distinct ordinal link functions/models, each with its own kernel and its own `smart_cold_start` branch; the current single "Ordinal Regression" row (OLS on `y`, quantile-mapped intercepts) only documents the proportional-odds kernel (`EDI/src/fast_ordinal_regression.cpp:167`) and should not be assumed to generalize. Sources: `EDI/src/fast_ordinal_glmm.cpp:271`, `EDI/src/fast_ordinal_cauchit_regression.cpp`, `EDI/src/fast_ordinal_cloglog_regression.cpp`, `EDI/src/fast_ordinal_probit_regression.cpp`, `EDI/src/fast_adjacent_category_logit.cpp:216`, `EDI/src/fast_continuation_ratio_regression.cpp:149`, `EDI/src/fast_stereotype_logit.cpp`.
- [ ] TODO-14: **GEE paths** (`InferenceCountKKGEE`, ordinal/incidence KK-GEE families) — add a table row noting that GEE's initial fit reuses the Logistic/Poisson smart-start heuristic unconditionally: `EDI/src/fast_gee.cpp:131-132` calls `fast_logistic_regression_internal`/`fast_poisson_regression_internal` with `smart_cold_start` hardcoded to `true` (not policy-gated the way the standalone Logistic/Poisson paths are per TODO-1 — GEE does not expose its own toggle at all).
