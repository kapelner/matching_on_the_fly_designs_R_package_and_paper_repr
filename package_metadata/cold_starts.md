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
| **Weibull AFT** | **OLS on $\ln(y) + 0.572$**: (R Standard) Uses uncensored observations with a Gumbel mean shift. | **Moment-based scale**: $\ln(\sigma) = 0.5 \ln(\text{resid\_var} / 1.64)$. |
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
