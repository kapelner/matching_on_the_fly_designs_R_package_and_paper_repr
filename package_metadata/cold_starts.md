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

## Benchmark: Naive (Zero) vs. Smart Cold Starts

The following benchmark compares the performance of "Naive" initialization (all parameters set to zero) against the "Smart Cold Start" strategies described above. 

**Benchmark Environment:**
* **Sample Size ($n$):** 1000
* **Covariates ($p$):** 10 (Linear combination, $\sum \beta^2 = 1$)
* **Replications:** 10 (using `microbenchmark`)
* **Hardware:** Linux (Ubuntu 15.2.0), Multi-core.

### Performance Summary

| Model | Naive It. | Smart It. | It. Improv. | Naive Time (ms) | Smart Time (ms) | Time Improv. |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Logistic** | 5 | 5 | 0.0% | 6.66 | 10.71 | -60.9% |
| **Poisson** | 6 | 6 | 0.0% | 4.44 | 9.95 | -123.9% |
| **Weibull** | 6 | 12 | -100.0% | 18.37 | 90.20 | -391.1% |
| **NegBin** | 19 | 4 | 78.9% | 18.12 | 12.85 | 29.1% |
| **Beta** | N/A | N/A | N/A | 170.13 | 122.01 | 28.3% |
| **Ordinal** | 7 | 6 | 14.3% | 601.57 | 568.59 | 5.5% |
| **Stereotype**| N/A | N/A | N/A | 3110.0 | 3460.0 | -11.4% |
| **ZINB** | N/A | N/A | N/A | 41.39 | 46.31 | -11.9% |
| **LogisticGLMM**| N/A | N/A | N/A | 833.24 | 1350.96 | -62.1% |
| **PoissonGLMM** | N/A | N/A | N/A | 193.87 | 446.06 | -130.1% |

### Analysis of Benchmark Results

The benchmark indicates that while "Smart Cold Starts" are transformative for certain model classes, they introduce performance overhead in others.

1.  **The Cost of "Smart" is Non-Zero**
    Computing a smart start (like an OLS fit) involves:
    *   Allocating memory for temporary matrices.
    *   Computing the cross-products $X^TX$ and $X^Ty$.
    *   Performing a matrix decomposition (QR or Cholesky).
    *   Transforming the response (log, logit, etc.).
    For a sample size of $n=1000$ and $p=10$, these operations take several milliseconds in C++.

2.  **"Well-Behaved" Data is Easy to Solve**
    In the benchmark, the data is generated from a latent linear model with standard normal covariates. This produces a "smooth" likelihood surface.
    *   **Logistic/Poisson:** As shown in the table, these models converge in only 5 to 6 iterations starting from zero.
    *   Because the naive start is already "close enough" to the global maximum, the time spent calculating the OLS solution is significantly greater than the time saved by reducing the iteration count by 1 or 2.

3.  **The "Double Fit" Penalty (Weibull)**
    The Weibull results are the most extreme (-391% slower) because of the intercept-only refinement.
    *   The Weibull "smart start" first does an OLS fit on $\ln(y)$, then performs a complete 1-parameter MLE fit for the intercept and scale to ensure the optimizer starts in a high-probability region.
    *   Running an entire sub-optimization before the main optimization is very slow for "easy" data, though it provides massive protection against divergence on "hard" data.

4.  **GLMM Overhead**
    In GLMMs (Logistic/Poisson), the smart start only initializes the fixed-effects $\beta$.
    *   The most expensive part of the GLMM likelihood is the Gauss-Hermite quadrature over the random effects.
    *   Doing an OLS fit for $\beta$ doesn't help the random effect variance parameter ($\sigma_u$) much, and the overhead of the OLS fit is amplified by the fact that the GLMM likelihood itself is already computationally heavy.

5.  **Why do we still use them?**
    If they are slower, why bother?
    *   **Stability:** On real-world datasets with separation (in logistic), extreme rates (in Poisson), or heavy censoring (in Weibull), a naive start of zero can often lead the optimizer to a region where the Hessian is singular or the gradient is flat, causing the model to fail to converge entirely.
    *   **The "Cold" vs "Warm" distinction:** These smart starts are only for the **first** subject in a sequential design (when you have 0 prior data). For the 2nd through 1000th subject, EDI uses **Warm Starts** (the result of the previous subject), which are nearly instantaneous and much faster than naive starts.

**Summary:** The benchmark shows that for a single "one-off" fit on easy data, naive starts are faster due to zero overhead. Smart starts are an insurance policy for robustness at the beginning of an experiment.

## Data Generation Models

The benchmark utilizes the synthetic data generation logic from `SimulationFramework.R`:

### 1. Latent Linear Model
Covariates $X$ are drawn from standard normal distributions. The latent continuous signal is generated as:
$$y_{cont} = X\beta$$
where $\beta$ is a vector of coefficients evenly spaced from 1 to -1, rescaled such that $\sum \beta^2 = 1$.

### 2. Response Transformations
The latent signal $y_{cont}$ is transformed based on the model family:
* **Incidence:** $\text{Bernoulli}(\text{plogis}(y_{cont}))$
* **Count:** $\text{Poisson}(\exp(y_{cont}))$ or $\text{NegBin}(\exp(y_{cont}), \theta=2)$
* **Survival:** $\text{Exp}(\exp(y_{cont}))$ with 25% independent censoring.
* **Proportion:** $\text{Beta}(\mu, \phi=10)$ where $\mu = \text{plogis}(y_{cont})$.
* **Ordinal:** Categorized into 4 levels based on quantiles of $y_{cont}$.

