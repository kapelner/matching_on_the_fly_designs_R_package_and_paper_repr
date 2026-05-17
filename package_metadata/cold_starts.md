# Smart Cold Start Strategies in EDI

This report summarizes how the "smart cold start" (initial parameter selection from no prior knowledge) is performed for each likelihood path in the `EDI` package. Strategies marked with "(R Standard)" have been aligned with core R functions like `glm.fit` and `survival::survreg`.

| Model / Likelihood Path | Primary Strategy for $\beta$ (Treatment & Covariates) | Strategy for Secondary Parameters ($\alpha, \sigma, \theta$, etc.) |
| :--- | :--- | :--- |
| **Logistic Regression** | **OLS on $\text{logit}(\frac{y + 0.5}{2})$**: (R Standard) Shrinks 0/1 response toward 0.5 to provide a robust initial separation plane. | N/A |
| **Poisson Regression** | **OLS on $\ln(y + 0.1)$**: (R Standard) Uses a small offset for zero counts to allow log-linear initialization. | N/A |
| **Negative Binomial** | **OLS on $\ln(y + 1)$**: Log-linear approximation using OLS. | $\ln(\theta) = 0$ (initial dispersion set to 1). |
| **Log-Binomial** | **OLS on $\ln(y + 1)$**: Or OLS on $y$ depending on the link function. | N/A |
| **Weibull AFT** | **OLS on $\ln(y) + 0.572$**: (R Standard) Uses uncensored observations with a Gumbel mean shift. | **Moment-based scale**: $\ln(\sigma) = 0.5 \ln(\text{resid\_var} / 1.64)$ (R Standard). |
| **Cox PH** | **OLS on $\ln(y)$**: Approximates log-hazard relative to an OLS baseline. | N/A |
| **Ordinal Regression** | **OLS on $y$**: Uses ordinal levels as a continuous response for $\beta$. | **Quantile Mapping**: $\alpha$ intercepts are set by mapping sample quantiles through the link's quantile function. |
| **Stereotype Logit** | **OLS on $y$**: Uses ordinal levels for initial $\beta$. | **Log-Frequency**: $\alpha$ intercepts based on category counts; category scores ($\phi$) start at 0 (equidistant). |
| **Beta Regression** | **OLS on $\text{logit}(y)$**: Maps (0,1) response to real line via logit. | $\ln(\phi) = \ln(\text{start\_phi})$ (typically 2.0). |
| **Zero-Inflated (ZINB)** | **Count**: OLS on $\ln(y+1)$. **ZI**: OLS on $I(y=0)$ indicator. | $\ln(\theta) = 0$. |
| **Zero-One Inflated Beta**| **Beta**: OLS on $\text{logit}(y)$ for $(0,1)$ subset. **Z/O**: OLS on $I(y=0)$ and $I(y=1)$. | $\ln(\phi) = 2.0$. |
| **GLMMs (Logistic/Poisson)** | **OLS on response or log-response**: Similar to fixed-effects versions. | **Small-Variance Start**: $\ln(\sigma_u)$ typically initialized to -3.0 or $\ln(\sigma_e/2)$. |
| Robust Regression      | **OLS on $y$**: Standard OLS to get an initial "clean" baseline. | Residual-based scale estimation. |

## 2026 Updated Benchmark: End-to-End Warm Start Flow

This updated benchmark evaluates the **full end-to-end process** required to use a smart start. Unlike theoretical benchmarks that compare iteration counts alone, this measures the total time including:
1.  **Heuristic Generation:** The time spent calculating OLS or moment-based starts (using `fast_ols_cpp`).
2.  **Argument Prep:** R-side overhead of repetitive vector construction and list passing.
3.  **Optimization:** The actual C++ solver execution time.

**Benchmark Environment:**
* **Sample Size ($n$):** 500
* **Covariates ($p$):** 5
* **Signal Strength:** Very Strong ($\beta \sim N(0, 2.0^2)$).
* **Replications:** 20 (using `system.time` for end-to-end flow).
* **Warm Start Heuristics:**
    *   **General:** Full warm start ($\beta$ from OLS, diagonal IRLS weights, and Hessian pre-calculation).
    *   **Logistic (IRLS):** Weight-only warm start derived from the mean success rate, mimicking `glm.fit`.

### Performance Summary (End-to-End)

| Model Family | Cold Start | Warm Start (Total) | Speedup (%) | Primary Benefit |
| :--- | :---: | :---: | :---: | :--- |
| **Negative Binomial** | 8.90 ms | 6.30 ms | **+29.2%** | OLS-based $\theta$ and $\beta$ jump-start. |
| **Logistic GLMM** | 36.15 ms | 25.70 ms | **+28.9%** | Reduced iterations in random-effect integration. |
| **Weibull (AFT)** | 0.45 ms | 0.40 ms | **+11.1%** | Improved initial scale ($\sigma$) estimate. |
| **Ordinal Regression** | 1.30 ms | 1.75 ms | *-34.6%* | Overhead of quantile mapping. |
| **Beta Regression** | 7.10 ms | 22.00 ms | *-209.9%* | Heuristic mismatch with L-BFGS surface. |
| **Logistic (IRLS)** | 0.15 ms | 2.50 ms | *-1567%* | R-side argument prep overhead. |

### Analysis of High-Signal Results

1.  **High-Complexity Wins**
    The most significant gains are found in **Negative Binomial** and **GLMM** models. For these "heavier" likelihoods, the $\sim 0.3\text{ms}$ cost of running `fast_ols_cpp` is negligible compared to the time saved during the iterative optimization. In **NegBin**, the OLS start provides a much better initial guess for the dispersion than a naive zero start.

2.  **The "Overhead Floor" for Simple Models**
    For models that are already extremely fast (converging in <1ms), the "Warm Start Flow" is actually a net negative due to R-side overhead.
    *   **Logistic (IRLS):** While the C++ solver is extremely efficient, the R-side overhead of calculating the mean and passing `warm_start_weights` adds significantly more time than the iterations saved.
    *   **Poisson (IRLS):** Shows a similar pattern where the cold start is so fast that even the fastest OLS routine adds unnecessary cost.

3.  **Robustness vs. Speed**
    Even in cases where the "Warm Start" is slower, it provides an **insurance policy**. A cold start at zero can fail to converge if the likelihood surface is non-concave or has regions of flat gradients. The Smart Cold Start ensures that every single inference path starts from a statistically motivated position, minimizing "Convergence Failure" errors returned to the user.

---

## Data Generation Models

The benchmark utilizes a "Very Strong Signal" synthetic data generation logic:

### 1. Latent Linear Model
Covariates $X$ are drawn from standard normal distributions ($p=5$). The latent continuous signal is generated as:
$$y_{cont} = X\beta$$
where $\beta \sim N(0, 2.0^2)$. This creates a highly varied likelihood surface that tests optimizer robustness.

### 2. Response Transformations
The latent signal $y_{cont}$ is transformed based on the model family:
* **Incidence:** $\text{Bernoulli}(\text{plogis}(y_{cont}))$
* **Count:** $\text{Poisson}(\exp(y_{capped}))$ or $\text{NegBin}(\exp(y_{capped}), \theta=2)$ (capped at $\pm 4$).
* **Survival:** $\text{Exp}(\exp(y_{capped}))$ with 20% independent censoring.
* **Proportion:** $\text{Beta}(\mu, \phi=10)$ where $\mu = \text{plogis}(y_{cont})$.
* **Ordinal:** Categorized into 4 levels based on proportional odds thresholds.
