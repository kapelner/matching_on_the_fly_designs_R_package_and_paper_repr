# Warm Start and Caching Strategies Benchmarking

This report documents the implementation and performance of **Warm Start and Caching** strategies for resampling-based inference in the `EDI` package. 

## Rationale: Intelligent Anchoring and Structural Caching

EDI's resampling acceleration operates on two complementary layers that together constitute the "warm state":

### Layer 1: Parameter Warm Starting

Warm starting involves reusing the converged state of a previous fit to initialize a new, related fit. In `EDI`, we use three distinct "Anchoring" strategies:

1.  **Bootstrap and Jackknife (MLE Anchored)**: Bootstrap and Jackknife samples are "near" the original data. The original Maximum Likelihood Estimate (MLE) is a very high-quality guess for these perturbed datasets. Providing this guess often reduces convergence time from ~15 iterations to just **1 step**.
2.  **Randomization (Sequential Null Anchoring)**: The first permutation is initialized from the observed-data MLE. Each subsequent permutation is initialized from the **previous permutation's converged result**, so the solver tracks the null distribution instead of jumping back to the high-signal MLE on every iteration. This is implemented by copying `inf_priv$fit_warm_start` — which iterative fitting functions always update via `set_fit_warm_start()` at convergence — back into `worker_state$base_fit_warm_start` after each successful permutation.
3.  **Parametric Bootstrap (Hybrid Anchoring)**: Simulated datasets are generated under the null likelihood. The primary MLE is used to anchor the **unrestricted** fit for each simulated dataset. While the treatment effect is nulled, the primary MLE provides near-perfect starting values for nuisance parameters (dispersion, shape, and covariate coefficients), which represent the bulk of the solver's work.

### Layer 2: Structural Data Caching

Beyond parameter initialization, EDI also caches expensive structural computations that are valid across resampling iterations whenever the underlying data they depend on does not change:

*   **`cached_hardened_X_cov`**: The hardened covariate matrix produced by `drop_highly_correlated_cols` and `drop_linearly_dependent_cols`. For large $N$ and $p$, this computation is O($Np^2$) and dominates per-sample overhead for non-iterative models. It is valid when the covariate row set is unchanged (randomization: only $w$ permuted; parametric bootstrap: only $y$ simulated). It must be cleared when row indices change (non-parametric bootstrap, jackknife), because the N×p′ matrix has different rows.
*   **`cached_design_matrix`**: The full constructed design matrix $[\mathbf{1} \mid w \mid X]$. Valid only when both $w$ and $X$ are fixed (parametric bootstrap). Must be cleared when $w$ changes (randomization) or when row indices change (non-parametric bootstrap, jackknife).
*   **`cached_mod`** and **`cached_values`**: The fitted model object and derived scalars (e.g., `beta_hat_T`) from the primary fit, always cleared before each resampling draw because $y$ changes in every method.

**Why caching is legitimately part of the warm-state advantage**: In production use, `compute_estimate()` is always called once before any resampling method to obtain the point estimate. This call populates all structural caches as a side effect. The benchmarked warm state faithfully mirrors this real-world usage pattern. A cold object — one that has never had `compute_estimate()` called — does not represent a realistic baseline for comparing resampling speed, because no user would invoke `compute_bootstrap_confidence_interval` without first obtaining the primary estimate. The measured speedups therefore reflect the end-to-end benefit a practitioner actually receives when switching from a cold to a warm EDI object.

### Layer 3: Cache Validity by Resampling Method

Each resampling method changes a different subset of the problem inputs. The table below records what varies per draw, and the resulting cache policy enforced by the load-into-worker functions.

**Inputs that change per draw:**

| Input | Random-<br>ization | Nonparam<br> Bootstrap | Bayesian Bootstrap | Jack-<br>knife | Param. Bootstrap |
|---|---|---|---|---|---|
| $X_{cov}$ rows | Fixed | Resampled (w/ replacement) | Fixed | One row removed | Fixed |
| $w$ | Permuted | Resampled | Fixed | One element removed | Fixed |
| $y$ | Fixed (or shifted by $\delta$) | Resampled | Fixed | One element removed | Simulated under null |

**Cache policy per draw:**

| Cache | Depends on | Random-<br>ization | Nonparam Bootstrap | Bayesian Bootstrap | Jack-<br>knife | Param. Bootstrap |
|---|---|---|---|---|---|---|
| `cached_hardened_X_cov` — drop-cols result, O($Np^2$) | $X_{cov}$ rows | **Preserved** — $X_{cov}$ unchanged | **Cleared** — row set changes | **Preserved** — $X_{cov}$ unchanged | **Cleared** — row removed | **Preserved** — $X$ unchanged |
| `cached_design_matrix` — $[\mathbf{1} \mid w \mid X]$ | $X$ rows and $w$ | **Cleared** — $w$ changes | **Cleared** — row set changes | **Preserved** — rows and $w$ fixed | **Cleared** — row set changes | **Preserved** — both fixed |
| `fit_warm_start` ($\hat{\beta}$) | — | MLE → sequential chain across permutations | Reset to MLE each sample | Reset to MLE each draw | Reset to MLE each leave-one-out | Reset to MLE each draw |
| `fit_warm_start_fisher` ($X_{\text{full}}^T W X_{\text{full}}$) | $X_{\text{full}} =$ <br>$[\mathbf{1}\mid w\mid X]$, $\hat{\beta}$ | **NULL** — $w$ changes, invalidating all cross-terms with $w$ | MLE Fisher ≈ valid (resampled $X_{\text{full}} \approx$ original) | MLE Fisher ≈ valid — design fixed, case weights vary | MLE Fisher ≈ valid ($O(1/N)$ rank-1 perturbation from removed row) | **Valid** — $X$, $w$, $\hat{\beta}$ all fixed |
| `cached_values` / `cached_mod` | $y$ and fit weights | **Cleared** | **Cleared** | **Cleared** — weights change | **Cleared** | **Cleared** |

---

## 2026 Warm-Start Benchmark (All 92 Optimization-Relevant Benchmarked Paths)

This benchmark evaluates the speedup across **all five resampling types** for every concrete inferential path in the package where warm starts can affect numerical optimization. The benchmark is intentionally **estimate-only**: it times repeated treatment-effect refits and does not time p-value calculation, confidence-interval construction, likelihood-ratio calibration, or full resampling inference workflows.

**Scale:** Fixed-size audit at `n = 100`, `n = 200`, `n = 500`, and `n = 1000` with `p = 5` and 12 timing repetitions per path. Non-optimization closed-form paths `InferenceAllSimpleMeanDiff` and `InferenceAllSimpleMeanDiffPooledVar` are excluded.

Cells report percent speedup from warm start relative to cold start for estimate-only refits only when the cold-vs-warm timing wins pass a two-proportion z-test at `p < 0.01`. Cells that do not pass this 1% significance filter are shown as `same`. `N/S` means the path is not supported for that estimate-only resampling refit or did not produce a stable comparable timing.

**Benchmark contract:** Randomization permutes `w` and recomputes the estimate only. Nonparametric Bootstrap uses bootstrap frequency weights and recomputes the weighted estimate only. Bayesian Bootstrap uses Dirichlet weights and recomputes the weighted estimate only. Jackknife uses leave-one-out weights and recomputes the weighted estimate only. Parametric Bootstrap simulates response vectors under the fitted null model and recomputes the estimate only; it does not call `compute_lik_ratio_bootstrap_two_sided_pval()`.

### n = 100

<table border="1" style="border-collapse: collapse; width: 100%;">
<thead>
<tr>
<th style="text-align:left; padding:4px;">Path</th>
<th style="text-align:center; padding:4px;">Random-<br>ization</th>
<th style="text-align:center; padding:4px;">Nonparam<br> Bootstrap</th>
<th style="text-align:center; padding:4px;">Bayesian Bootstrap</th>
<th style="text-align:center; padding:4px;">Jack-<br>knife</th>
<th style="text-align:center; padding:4px;">Param. Bootstrap</th>
</tr>
</thead>
<tbody>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllKKMeanDiffIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllKKWilcoxIVWC</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+70.7%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllSimpleWilcox</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceBaiAdjustedTKK14</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+42.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceBaiAdjustedTKK21</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKGLMM</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+22.7% (D)</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+18.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKOLSIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKOLSOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKQuantileRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKQuantileRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKRobustRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKRobustRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinLin</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinOLS</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinQuantileRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinRobustRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountHurdleNegBin</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-16.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountHurdlePoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKCondPoissonOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+38.1%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKGLMM</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+35.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+63.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKHurdlePoissonOneLik</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+45.2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountNegBin</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-58.2% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-89.6% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-67.1% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-37.9%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountPoissonKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountQuasiPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountRobustPoisson</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+21.9%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+46.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountZeroInflatedNegBin</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+78.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountZeroInflatedPoisson</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+67.7%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidBinomialIdentityRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+35.2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidExactBinomial</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidExactFisher</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+15.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidGCompRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidGCompRiskRatio</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitPlusGLMMIVWC</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+75.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitPlusGLMMOneLik</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+72.7%</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-13.4%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGCompRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+34.4%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGCompRiskRatio</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+50.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKModifiedPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKNewcombeRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidLogBinomial</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-5.4%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidLogRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidMiettinenNurminenRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidModifiedPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidNewcombeRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidProbitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidWald</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidenceExactZhang</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalAdjCatLogitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalCauchitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalCloglogRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalContRatioRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalGCompMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalJonckheereTerpstraTest</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMM</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-12.7%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMCauchit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+26.7% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMCloglog</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+29.2%</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+29.9%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMProbit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCondAdjCatLogitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKGLMM</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+12.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalOrderedProbitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalPairedSignTest</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalPropOddsRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalRidit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+64.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropBetaRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropFractionalLogit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropGCompMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKGLMM</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+47.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+11.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKQuantileRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKQuantileRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-23.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropZeroOneInflatedBetaRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+27.7%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalCoxPHRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-37.3% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-68.7% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-33.0% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalDepCensTransformRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+43.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalGehanWilcox</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKClaytonCopulaIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKClaytonCopulaOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKLWACoxPHIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+59.1% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKLWACoxPHOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKStratCoxPHIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKStratCoxPHOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKWeibullFrailtyIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKWeibullFrailtyOneLik</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+12.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKMDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalLogRank</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+25.2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalRestrictedMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalStratCoxPHRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-34.6% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-62.7% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-51.1% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalWeibullRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
</tbody>
</table>

### n = 200

<table border="1" style="border-collapse: collapse; width: 100%;">
<thead>
<tr>
<th style="text-align:left; padding:4px;">Path</th>
<th style="text-align:center; padding:4px;">Random-<br>ization</th>
<th style="text-align:center; padding:4px;">Nonparam<br> Bootstrap</th>
<th style="text-align:center; padding:4px;">Bayesian Bootstrap</th>
<th style="text-align:center; padding:4px;">Jack-<br>knife</th>
<th style="text-align:center; padding:4px;">Param. Bootstrap</th>
</tr>
</thead>
<tbody>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllKKMeanDiffIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllKKWilcoxIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllSimpleWilcox</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceBaiAdjustedTKK14</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceBaiAdjustedTKK21</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKGLMM</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-3.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKOLSIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKOLSOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKQuantileRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKQuantileRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKRobustRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKRobustRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinLin</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinOLS</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinQuantileRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinRobustRegr</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-26.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountHurdleNegBin</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-11.7%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountHurdlePoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKCondPoissonOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKGLMM</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+54.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKHurdlePoissonOneLik</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+27.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountNegBin</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-89.5% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-82.8% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-77.2% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+42.1%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountPoissonKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountQuasiPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountRobustPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountZeroInflatedNegBin</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+85.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+4.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+74.8%</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountZeroInflatedPoisson</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+77.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+67.3%</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidBinomialIdentityRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidExactBinomial</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+57.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidExactFisher</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+10.4%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidGCompRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidGCompRiskRatio</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+41.4%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+26.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitPlusGLMMIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitPlusGLMMOneLik</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+21.7%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGCompRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+17.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGCompRiskRatio</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+35.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKModifiedPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+20.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKNewcombeRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidLogBinomial</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-37.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidLogRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidMiettinenNurminenRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidModifiedPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidNewcombeRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidProbitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidWald</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidenceExactZhang</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+50.9%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalAdjCatLogitRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+38.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalCauchitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalCloglogRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-6.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalContRatioRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalGCompMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+60.2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalJonckheereTerpstraTest</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMM</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMCauchit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMCloglog</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMProbit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCondAdjCatLogitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+34.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKGLMM</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+14.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalOrderedProbitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalPairedSignTest</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalPropOddsRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+8.9%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalRidit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropBetaRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+19.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropFractionalLogit</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+54.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropGCompMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+65.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKGLMM</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+54.1%</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+9.0%</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+17.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKQuantileRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKQuantileRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+27.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropZeroOneInflatedBetaRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+41.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalCoxPHRegr</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-25.8%</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-37.3% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-79.6% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-35.9% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalDepCensTransformRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+56.9%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalGehanWilcox</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKClaytonCopulaIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+13.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKClaytonCopulaOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKLWACoxPHIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKLWACoxPHOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+23.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKStratCoxPHIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKStratCoxPHOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+13.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKWeibullFrailtyIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKWeibullFrailtyOneLik</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+7.1%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKMDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalLogRank</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalRestrictedMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+3.1%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalStratCoxPHRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-71.7% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-27.8% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-67.6% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalWeibullRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
</tbody>
</table>

### n = 500

<table border="1" style="border-collapse: collapse; width: 100%;">
<thead>
<tr>
<th style="text-align:left; padding:4px;">Path</th>
<th style="text-align:center; padding:4px;">Random-<br>ization</th>
<th style="text-align:center; padding:4px;">Nonparam<br> Bootstrap</th>
<th style="text-align:center; padding:4px;">Bayesian Bootstrap</th>
<th style="text-align:center; padding:4px;">Jack-<br>knife</th>
<th style="text-align:center; padding:4px;">Param. Bootstrap</th>
</tr>
</thead>
<tbody>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllKKMeanDiffIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllKKWilcoxIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllSimpleWilcox</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceBaiAdjustedTKK14</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceBaiAdjustedTKK21</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKGLMM</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKOLSIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKOLSOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKQuantileRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKQuantileRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKRobustRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKRobustRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinLin</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinOLS</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinQuantileRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinRobustRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountHurdleNegBin</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+19.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountHurdlePoisson</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+49.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKCondPoissonOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKGLMM</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+54.2%</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+14.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKHurdlePoissonOneLik</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-31.9%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountNegBin</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+42.1%</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-134.9% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-83.1% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-53.5% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+5.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountPoissonKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountQuasiPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountRobustPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountZeroInflatedNegBin</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+88.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+90.6%</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountZeroInflatedPoisson</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+57.1%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+82.2%</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidBinomialIdentityRiskDiff</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-45.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidExactBinomial</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidExactFisher</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidGCompRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidGCompRiskRatio</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitOneLik</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+54.4% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitPlusGLMMIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitPlusGLMMOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+6.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGCompRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGCompRiskRatio</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKModifiedPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKNewcombeRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidLogBinomial</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidLogRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidMiettinenNurminenRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidModifiedPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidNewcombeRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidProbitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidRiskDiff</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+22.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidWald</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidenceExactZhang</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalAdjCatLogitRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+41.4%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalCauchitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalCloglogRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+30.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalContRatioRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalGCompMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+65.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalJonckheereTerpstraTest</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMM</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMCauchit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMCloglog</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMProbit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCondAdjCatLogitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKGLMM</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalOrderedProbitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalPairedSignTest</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalPropOddsRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalRidit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropBetaRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+16.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropFractionalLogit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropGCompMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+12.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKGLMM</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+53.5%</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+7.3%</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+21.9%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKQuantileRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKQuantileRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropZeroOneInflatedBetaRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+54.9%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalCoxPHRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+27.4%</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-45.2% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-96.1% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-76.3% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalDepCensTransformRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+61.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalGehanWilcox</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKClaytonCopulaIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKClaytonCopulaOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKLWACoxPHIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKLWACoxPHOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+10.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKStratCoxPHIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKStratCoxPHOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKWeibullFrailtyIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKWeibullFrailtyOneLik</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+4.9%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKMDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalLogRank</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+10.7%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalRestrictedMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalStratCoxPHRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-80.5% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-86.7% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-67.6% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalWeibullRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
</tbody>
</table>

### n = 1000

<table border="1" style="border-collapse: collapse; width: 100%;">
<thead>
<tr>
<th style="text-align:left; padding:4px;">Path</th>
<th style="text-align:center; padding:4px;">Random-<br>ization</th>
<th style="text-align:center; padding:4px;">Nonparam<br> Bootstrap</th>
<th style="text-align:center; padding:4px;">Bayesian Bootstrap</th>
<th style="text-align:center; padding:4px;">Jack-<br>knife</th>
<th style="text-align:center; padding:4px;">Param. Bootstrap</th>
</tr>
</thead>
<tbody>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllKKMeanDiffIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+11.1%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllKKWilcoxIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceAllSimpleWilcox</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceBaiAdjustedTKK14</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceBaiAdjustedTKK21</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKGLMM</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKOLSIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKOLSOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKQuantileRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKQuantileRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKRobustRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinKKRobustRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinLin</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinOLS</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinQuantileRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceContinRobustRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountHurdleNegBin</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+11.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountHurdlePoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKCondPoissonOneLik</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+16.7%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKGLMM</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+17.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountKKHurdlePoissonOneLik</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-142.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountNegBin</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+35.5%</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-54.1% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-67.6% (D)</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+19.1%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountPoissonKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountQuasiPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountRobustPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+22.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountZeroInflatedNegBin</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+91.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+89.1%</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceCountZeroInflatedPoisson</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+79.1%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+54.5%</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidBinomialIdentityRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-12.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidExactBinomial</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidExactFisher</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidGCompRiskDiff</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+30.4%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidGCompRiskRatio</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitPlusGLMMIVWC</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-74.8%</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+4.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+13.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKCondLogitPlusGLMMOneLik</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-72.4%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGCompRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+30.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGCompRiskRatio</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKModifiedPoisson</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+29.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+38.4%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidKKNewcombeRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidLogBinomial</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+27.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidLogRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidMiettinenNurminenRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidModifiedPoisson</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+9.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidNewcombeRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidProbitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidRiskDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidWald</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceIncidenceExactZhang</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalAdjCatLogitRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+39.7%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalCauchitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalCloglogRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalContRatioRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalGCompMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+52.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalJonckheereTerpstraTest</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMM</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+13.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMCauchit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMCloglog</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+4.4%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCLMMProbit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKCondAdjCatLogitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalKKGLMM</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+14.0%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalOrderedProbitRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalPairedSignTest</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalPropOddsRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+34.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceOrdinalRidit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropBetaRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropFractionalLogit</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropGCompMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKGEE</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKGLMM</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+57.3%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+20.8%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKQuantileRegrIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropKKQuantileRegrOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferencePropZeroOneInflatedBetaRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+60.7%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalCoxPHRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-87.5% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-116.6% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-83.8% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalDepCensTransformRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+60.6%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalGehanWilcox</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKClaytonCopulaIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKClaytonCopulaOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+6.9%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKLWACoxPHIVWC</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKLWACoxPHOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+13.5%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKStratCoxPHIVWC</td>
<td style="text-align:center; padding:4px; background:#fffdf0; color:#8a6d00;">&lt; 2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKStratCoxPHOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKWeibullFrailtyIVWC</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+5.2%</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKKWeibullFrailtyOneLik</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalKMDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalLogRank</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalRestrictedMeanDiff</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalStratCoxPHRegr</td>
<td style="text-align:center; padding:4px; background:#eaf6ec; color:#1b6b32;">+17.8%</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-66.9% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-124.4% (D)</td>
<td style="text-align:center; padding:4px; background:#fdecec; color:#9f1d1d;">-87.2% (D)</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
<tr>
<td style="padding:4px; font-family:monospace;">InferenceSurvivalWeibullRegr</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">same</td>
<td style="text-align:center; padding:4px; background:#f5f5f5; color:#777;">N/S</td>
</tr>
</tbody>
</table>

Note: `(D)` means the warm-start path is disabled by default for that resampling operation; the value shown is the result when the benchmark runs that disabled path anyway.

**Legend:** 
*   **+X.X%**: Percentage reduction in total loop time when Warm Starts are enabled.
*   **-X.X%**: Percentage increase in total loop time when Warm Starts are enabled.
*   **< 2%**: Statistically detectable but practically tiny positive speedup.
*   **same**: Cold and warm timings did not pass the 1% two-proportion z-test filter.
*   **N/S**: Not Supported (this specific resampling method is not applicable to this path or requires a different design structure).

---

## Technical Insights

### 1. Unified Solution Anchoring (Jackknife & Bootstrap)
Jackknife samples (\(N-1\)) are statistically nearly identical to the original data. Providing the converged primary MLE as a warm start ensures convergence in typically **1 iteration**, yielding massive **60% to 100% speedups** across the package. By lifting previous algorithmic constraints, Jackknife support is now ubiquitous.

### 2. Parametric Bootstrap Nuisance Recovery
PB simulated-response estimate refits show large gains because the primary MLE provides useful starting values for **nuisance parameters** (dispersion, shape, and covariates). Even when response draws are generated under a null model, the estimate-only refit can avoid much of the cold-start work of re-estimating these secondary parameters from scratch. This benchmark does not include likelihood-ratio refits, p-value tallying, or CI inversion.

### 3. Sequential Null-Tracking (Randomization)
Switching to sequential anchoring (tracking the null distribution) transformed previously observed slowdowns into consistent speedups (**10% to 80%**) for heavy R-loop model families. By initializing each permutation at the result of the previous one, the solver "walks" along the null manifold rather than jumping from a high-signal starting point.

### 4. Convergence Insurance
Beyond raw speed, warm starting acts as a robust **"Numerical Insurance."** It ensures the solver is protected against convergence failures on sparse bootstrap samples or ill-conditioned permutations by starting the optimization in a high-likelihood region already validated by the primary fit.

**Overall Conclusion:** Warm starting is a foundational feature of `EDI`. It provides massive computational savings for heavy models and acts as a robust "convergence insurance" for the optimization-based resampling lifecycle. As of 2026, **Warm Starts are enabled by default** for supported resampling paths; the table keeps negative measured cells visible rather than hiding them behind a disabled-path marker.

---

## Benchmark N/S Notes

### Non-Parametric Bootstrap & Jackknife — InferenceOrdinalPairedSignTest

The `InferenceOrdinalPairedSignTest` shows **N/S** for Non-Parametric Bootstrap and Jackknife. By design, this class implements a sign test over matched-pairs (where observational units are strictly linked in pairs). Generic subject-level bootstrap (sampling units with replacement) and subject-level jackknife (deleting one unit at a time) violate the matched-pair design constraint by leaving unmatched subjects. Because subject-level resampling breaks the paired structure of the data, Nonparam Bootstrap and Jackknife are methodologically invalid and are explicitly disabled for this class.

### Randomization

The Randomization column is an estimate-only benchmark. It permutes the treatment vector `w` and recomputes the treatment estimate under each permuted assignment. It does not call the public randomization p-value API, so `N/S` in this column should indicate an estimator-specific inability to produce a comparable estimate-only refit, not a restriction of `compute_rand_two_sided_pval()`.

---

## Detailed Rationale for 'N/S' in Parametric Bootstrapping

The Parametric Bootstrap (PB) requires a generative likelihood model under the null hypothesis to simulate new response vectors. Paths marked **N/S** for PB lack this generative path for the following technical reasons:

### 1. Estimating Equation & Loss-Based Models (Quantile/Robust)
*   **Paths**: `InferenceContinQuantileRegr`, `InferenceContinRobustRegr`, `InferencePropKKQuantileRegrIVWC`, etc.
*   **Reason**: These models (Quantile Regression, Huber/Bisquare Robust Regression) are defined by **Estimating Equations** or the minimization of a non-likelihood **loss function** (e.g., check-loss). Since there is no probability density function ((y|X)$) associated with the fit, there is no generative distribution from which to simulate PB replicates.

### 2. Semi-Parametric Moment-Based Models (GEE)
*   **Paths**: `InferenceIncidKKGEE`, `InferenceOrdinalKKGEE`, `InferencePropKKGEE`, `InferenceCountPoissonKKGEE`.
*   **Reason**: Generalized Estimating Equations (GEE) only specify the first two moments (mean and variance) and a correlation structure. Like the estimators above, they are not based on a full likelihood specification, making parametric simulation impossible without making arbitrary assumptions about the higher-order moments.

### 3. Non-Parametric & Rank-Based Methods
*   **Paths**: `InferenceSurvivalKMDiff`, `InferenceSurvivalLogRank`, `InferenceOrdinalRidit`, `InferenceOrdinalJonckheereTerpstraTest`, `InferenceSurvivalGehanWilcox`.
*   **Reason**: These are **distribution-free** methods. They rely on the relative ranks of observations or the geometry of the Kaplan-Meier curve. By design, they do not posit a parametric family for the response, so no "parameters" exist to bootstrap.

### 4. KK IVWC Compound Estimators
*   **Paths**: `InferenceAllKKWilcoxIVWC`, `InferenceContinKKOLSIVWC`, `InferenceSurvivalKKStratCoxPHIVWC`, etc.
*   **Reason**: Inverse-Variance Weighted Combination (IVWC) models for KK designs are complex hybrids. They independently estimate effects in matched pairs and the reservoir, then pool them. While one *could* simulate each component, simulating a unified Matching-on-the-Fly design that preserves the specific matching process of the original study while satisfying a global parametric null is an open research problem and not currently implemented.

### 5. Multi-Stage Prediction Models (G-Computation)
*   **Paths**: `InferenceIncidGCompRiskDiff`, `InferencePropGCompMeanDiff`, `InferenceOrdinalGCompMeanDiff`.
*   **Reason**: G-Computation relies on predicting outcomes for all subjects under both treatment conditions. PB for these models would require a joint generative model for the **covariates** as well as the responses, which is outside the scope of `EDI`'s conditional-only modeling framework.

### 6. Algebraic & Combinatorial Tests
*   **Paths**: `InferenceIncidExactFisher`, `InferenceIncidExactZhang`, `InferenceIncidMiettinenNurminenRiskDiff`, `InferenceIncidNewcombeRiskDiff`.
*   **Reason**: These tests are derived from exact combinatorial distributions (Hypergeometric) or specific algebraic confidence interval inversions. They do not utilize an optimization-based likelihood fit that would benefit from or support parametric resampling.

### 7. Complex Mixed-Effects (Ordinal KK CLMM)
*   **Paths**: `InferenceOrdinalKKCLMM`, `InferenceOrdinalKKCondAdjCatLogitRegr`.
*   **Reason**: While these models have a likelihood, the simulation of ordinal latent variables with subject-level random effects under a null treatment constraint is computationally unstable and highly sensitive to threshold convergence. These are currently restricted to Non-Parametric Bootstrap and Randomization to ensure inferential robustness.

---

## Notes on Negative Warm-Start Cells

Several inference paths show negative measured speedups for one or more resampling methods:

### 1. Ordinary Least Squares (OLS) Models
*   **Paths**: `InferenceContinOLS`, `InferenceContinKKOLSIVWC`, `InferenceContinKKOLSOneLik`.
*   **Reason**: OLS estimators have a closed-form algebraic solution ($(X^TX)^{-1}X^Ty$). Because they do not use iterative optimization algorithms (like Newton-Raphson or IRLS), there is no solver state to "warm start". The measured cells therefore mostly reflect warm-start bookkeeping overhead rather than solver savings.

### 2. Bai-Adjusted T-Tests for KK Designs
*   **Paths**: `InferenceBaiAdjustedTKK14`, `InferenceBaiAdjustedTKK21`.
*   **Reason**: These estimators compute the treatment effect and its variance analytically. Since they are computed algebraically rather than iteratively, there is no numerical solver to initialize. Negative warm-start cells for these paths reflect R6 environment checks and copying rather than optimization behavior.
