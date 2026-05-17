# Warm Starts in the EDI Package

This report documents the implementation and rationale for warm start strategies across different inferential procedures in the EDI package. Warm starting involves capturing the full converged state of a numerical optimizer (Parameters, IRLS Weights, and Fisher Information) and injecting it into subsequent, related fits to reduce iteration counts and total execution time.

## 1. Randomization Sampling (Permutation Tests)
*   **The Anchor**: The package performs one "observed fit" on the real data to find the Maximum Likelihood Estimate (MLE).
*   **The Strategy**: Every permutation (e.g., all 1,000 iterations) starts with the **Observed MLE** as its warm start.
*   **Why it works**: In a randomization test, only the treatment assignments are shuffled; the covariates and the outcome values remain the same. The "true" coefficients for the covariates from the observed fit are usually extremely close to the values needed for the permuted fits. This typically reduces the number of iterations per permutation from ~10 down to **0 or 1**.

## 2. Nonparametric Bootstrap Sampling
*   **The Anchor**: The **Observed MLE**.
*   **The Strategy**: Each bootstrap resample (drawn with replacement) is initialized using the coefficients and Hessian from the full dataset.
*   **Why it works**: A bootstrap resample is a "representative" subset of the original data. While some observations are doubled and others omitted, the overall likelihood surface is very similar. Starting at the full-data MLE puts the solver immediately in the "neighborhood" of the bootstrap MLE, avoiding the expensive initial search.

## 3. Parametric Bootstrap Sampling
*   **The Anchor**: The **Null Model** parameters.
*   **The Strategy**: In parametric bootstrapping (used for Likelihood Ratio Test calibration), we simulate new outcomes ($y^*$) based on the model fit under the null hypothesis ($H_0$). 
*   **Application**: When we fit the "Full" model on this simulated data, we warm-start it using the **Observed Null Parameters**.
*   **Why it works**: Since the data was generated from the null model, the simulated MLEs will be naturally clustered around the null values.

## 4. Jackknife Sampling
*   **The Anchor**: The **Observed MLE**.
*   **The Strategy**: The Jackknife performs $N$ fits, each omitting exactly one observation (or cluster). 
*   **Application**: Each "leave-one-out" fit is warm-started using the results of the $N$-sample fit.
*   **Why it works**: Omitting a single data point changes the likelihood surface by only $1/N$. The solution for $N-1$ points is mathematically almost identical to the solution for $N$ points. These fits often converge in **zero iterations** (the first gradient check passes the tolerance immediately).

## 5. Confidence Interval Inversion (Likelihood Profiles)
*   **The Anchor**: The **Previous Step** parameters.
*   **The Strategy**: Finding a CI via inversion involves a "root-finding" search (like bisection) across different values of the treatment effect ($\delta$).
*   **Application**: 
    1. Step 1: Fit $\delta = 0$. Cache result.
    2. Step 2: Fit $\delta = 0.1$. Warm-start with results from Step 1.
    3. Step 3: Fit $\delta = 0.15$. Warm-start with results from Step 2.
*   **Why it works**: This is a "sequential" warm start. As the root-finder narrows in on the confidence boundary, each guess is only a small perturbation from the previous one, allowing the solver to "track" the moving MLE with minimal effort.

---

## Summary Table: Warm Start Usage by Procedure

| Procedure | Initialization Anchor | Iteration Reduction |
| :--- | :--- | :--- |
| **Randomization** | Full Observed MLE | **High**: (Covariates are fixed) |
| **Jackknife** | Full Observed MLE | **Extreme**: (Data change is minimal) |
| **Nonparametric Boot** | Full Observed MLE | **Medium**: (Data change is moderate) |
| **Parametric Boot** | Null Model Parameters | **Medium**: (New data, but known source) |
| **CI Inversion** | Result of previous search step | **High**: (Sequential refinement) |

By implementing the **Tiered Policy**, the package ensures that "Heavy" models (like GLMMs) pass the full Hessian/Hessian bundle for these procedures, while "Light" models only pass the Beta vector to keep the R/C++ overhead low.
