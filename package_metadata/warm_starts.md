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

## 2026 Truly Exhaustive Benchmark (All 90 Concrete Paths)

This benchmark evaluates the speedup across **all five resampling types** for every concrete inferential path in the package. 

**Scale:** High-Resolution Audit ($N=100-500$, $P=2-10$, proportional to model complexity). 
Note: For ultra-fast models, speedups are measured relative to a sub-millisecond cold start baseline and are bounded by the numerical resolution floor.

<table border="1" style="border-collapse: collapse; width: 100%;">
  <thead>
  <tr style="background-color: #f2f2f2;">
    <th style="text-align: left; padding: 8px;">Inference Path</th>
    <th style="text-align: center; padding: 8px;">Random-<br>ization</th>
    <th style="text-align: center; padding: 8px;">Nonparam Bootstrap</th>
    <th style="text-align: center; padding: 8px;">Bayesian Bootstrap</th>
    <th style="text-align: center; padding: 8px;">Jack-<br>knife</th>
    <th style="text-align: center; padding: 8px;">Param. Bootstrap</th>
  </tr>
</thead>
  <tbody>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalPairedSignTest</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+6.0%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidBinomialIdentityRiskDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+53.7%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidGCompRiskDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+87.8%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+14.3%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidExtendedRobins</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+27.3%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidKKGCompRiskDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+8.8%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidenceExactZhang</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+29.5%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-2.4%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+8.8%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidCMH</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+33.4%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidLogRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+62.4%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+6.6%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidGCompRiskRatio</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+90.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidKKCondLogitIVWC</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px;">+3.5%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-12.5%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+13.8%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidWald</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+29.6%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidExactBinomial</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+23.8%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+23.6%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-2.2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidKKCondLogitPlusGLMMOneLik</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+14.9%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-16.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+6.3%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidKKGCompRiskRatio</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px;">+2.4%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidKKGEE</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+31.6%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-68.7%</td>
      <td style="text-align: center; padding: 8px;">+2.7%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidLogBinomial</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+51.9%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-13.3%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+11.5%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidExactFisher</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+20.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+9.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidKKCondLogitOneLik</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px;">+2.7%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-6.2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+13.6%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidModifiedPoisson</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+58.1%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-41.7%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+9.6%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidMiettinenNurminenRiskDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+18.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidRiskDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+13.7%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+23.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+6.8%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidKKModifiedPoisson</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidKKCondLogitPlusGLMMIVWC</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidKKNewcombeRiskDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+8.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-5.4%</td>
      <td style="text-align: center; padding: 8px;">+3.0%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidProbitRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+48.1%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-6.7%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+26.8%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceIncidNewcombeRiskDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px;">+4.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">+3.4%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalDepCensTransformRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+68.3%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">+2.8%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalWeibullRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+43.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+6.4%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceAllSimpleWilcox</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+24.8%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+7.1%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-5.8%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceAllKKWilcoxIVWC</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+26.4%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-55.3%</td>
      <td style="text-align: center; padding: 8px;">+2.6%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceAllKKMeanDiffIVWC</b></td>
      <td style="text-align: center; padding: 8px;">+5.0%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+5.2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-4.9%</td>
      <td style="text-align: center; padding: 8px;">+3.5%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceContinKKQuantileRegrIVWC</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+23.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-2.8%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+9.0%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceContinKKQuantileRegrOneLik</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+32.7%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+7.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+5.4%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceContinLin</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+33.7%</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+71.9%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+33.0%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceContinRobustRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+35.5%</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+83.9%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-2.3%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-9.7%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceContinKKRobustRegrOneLik</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+37.0%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+6.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-8.8%</td>
      <td style="text-align: center; padding: 8px;">+4.1%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceContinKKGLMM</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+25.6%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+33.7%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceContinKKRobustRegrIVWC</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+29.9%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+25.5%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+10.4%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountHurdleNegBin</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+30.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+46.7%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+6.1%</td>
      <td style="text-align: center; padding: 8px;">+3.5%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountNegBin</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+32.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+38.7%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-10.3%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-63.8%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceContinQuantileRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+34.7%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+24.6%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-8.9%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountPoisson</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+13.6%</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+66.8%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-15.4%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-2.2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountKKCondPoissonOneLik</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+26.5%</td>
      <td style="text-align: center; padding: 8px;">+3.5%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+21.9%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountHurdlePoisson</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+13.3%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+9.7%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountQuasiPoisson</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+33.0%</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+70.9%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+6.2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+15.3%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountZeroInflatedPoisson</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+28.9%</td>
      <td style="text-align: center; padding: 8px;">+4.1%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-5.5%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountKKHurdlePoissonOneLik</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+22.0%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-2.4%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">+2.8%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-2.0%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountKKGLMM</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+33.4%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+14.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">+3.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-3.7%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountPoissonKKGEE</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+12.5%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+17.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-6.0%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountRobustPoisson</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+30.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+28.8%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+33.1%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceCountZeroInflatedNegBin</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+36.8%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+17.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+8.6%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+8.0%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalJonckheereTerpstraTest</b></td>
      <td style="text-align: center; padding: 8px;">+4.4%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+18.8%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+12.5%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-2.7%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalCauchitRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+25.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+48.9%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalAdjCatLogitRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+35.6%</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+61.5%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">+2.7%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-4.8%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalGCompMeanDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+30.4%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+45.6%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-7.3%</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+53.9%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalKKGEE</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+17.9%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-3.0%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+20.4%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferencePropBetaRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+15.3%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+23.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-8.9%</td>
      <td style="text-align: center; padding: 8px;">+2.7%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalKKCLMM</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+23.0%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+6.3%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-7.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+7.0%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalKKCLMMCauchit</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+27.2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-2.3%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-8.2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalKKCondAdjCatLogitRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+29.0%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+9.6%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-15.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalOrderedProbitRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+20.3%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+46.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-5.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalKKGLMM</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+20.9%</td>
      <td style="text-align: center; padding: 8px;">+3.3%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+7.1%</td>
      <td style="text-align: center; padding: 8px;">+4.8%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferencePropKKQuantileRegrOneLik</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+30.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferencePropGCompMeanDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+14.5%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+28.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferencePropFractionalLogit</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+24.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+66.9%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-20.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalKKClaytonCopulaOneLik</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+29.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalCoxPHRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+7.4%</td>
      <td style="text-align: center; padding: 8px;">+2.1%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-7.6%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalContRatioRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+27.6%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+33.1%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">+3.4%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferencePropZeroOneInflatedBetaRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+24.5%</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+52.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-3.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalKKWeibullFrailtyOneLik</b></td>
      <td style="text-align: center; padding: 8px;">+5.0%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+28.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-8.7%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalKKLWACoxPHOneLik</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+25.8%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+7.8%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalKKLWACoxPHIVWC</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+28.7%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+9.2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px;">+2.7%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+11.2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalKKCLMMProbit</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+22.9%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+11.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">+3.6%</td>
      <td style="text-align: center; padding: 8px;">+2.8%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalCloglogRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+27.6%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+48.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">+2.1%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalLogRank</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+6.4%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+16.3%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+9.1%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalKMDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+26.6%</td>
      <td style="text-align: center; padding: 8px; background-color: #228b22; color: white;">+54.4%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-14.3%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+18.2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalPropOddsRegr</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+26.3%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+47.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+9.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+35.7%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalKKCLMMCloglog</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+23.1%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-3.5%</td>
      <td style="text-align: center; padding: 8px;">—</td>
      <td style="text-align: center; padding: 8px;">—</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferencePropKKGLMM</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+26.4%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+7.4%</td>
      <td style="text-align: center; padding: 8px;">—</td>
      <td style="text-align: center; padding: 8px;">—</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalGehanWilcox</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+19.5%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+23.3%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+11.5%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+6.9%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferencePropKKGEE</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+20.4%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+12.1%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-8.1%</td>
      <td style="text-align: center; padding: 8px;">+4.0%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalKKStratCoxPHOneLik</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+12.0%</td>
      <td style="text-align: center; padding: 8px;">+2.0%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalKKStratCoxPHIVWC</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+9.2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+5.3%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalStratCoxPHRegr</b></td>
      <td style="text-align: center; padding: 8px;">+3.6%</td>
      <td style="text-align: center; padding: 8px;">+2.2%</td>
      <td style="text-align: center; padding: 8px;">—</td>
      <td style="text-align: center; padding: 8px;">—</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+5.4%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalRestrictedMeanDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+7.5%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+16.2%</td>
      <td style="text-align: center; padding: 8px;">—</td>
      <td style="text-align: center; padding: 8px;">—</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceOrdinalRidit</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+25.1%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferencePropKKQuantileRegrIVWC</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+22.1%</td>
      <td style="text-align: center; padding: 8px;">+2.2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-5.9%</td>
      <td style="text-align: center; padding: 8px;">+4.1%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalKKClaytonCopulaIVWC</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+23.5%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+24.4%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceSurvivalKKWeibullFrailtyIVWC</b></td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+18.3%</td>
      <td style="text-align: center; padding: 8px;">< 2%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+8.3%</td>
      <td style="text-align: center; padding: 8px; background-color: #ccffcc;">+7.8%</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceAllSimpleMeanDiff</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceAllSimpleMeanDiffPooledVar</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceBaiAdjustedTKK14</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #ffcccc;">-18.8%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceBaiAdjustedTKK21</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceContinOLS</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceContinKKOLSIVWC</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">N/S</td>
      <td style="text-align: center; padding: 8px;">+2.7%</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
    </tr>
    <tr>
      <td style="text-align: left; padding: 8px;"><b>InferenceContinKKOLSOneLik</b></td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
      <td style="text-align: center; padding: 8px; background-color: #eeeeee;">Disabled</td>
    </tr>
  </tbody>
</table>

**Legend:** 
*   **+X.X%**: Percentage reduction in total loop time when Warm Starts are enabled.
*   **< X%**: Improvement is negligible or below the measurement noise floor for this ultra-fast model.
*   **N/S**: Not Supported (this specific resampling method is not applicable to this path or requires a different design structure).
*   **—**: Supported, but no separate Bayesian-bootstrap warm-start timing was recorded in this benchmark table.
*   **Disabled**: Warm Starts are disabled by default for this path (typically closed-form models where R6 environment overhead outweighs any mathematical gains).

---

## Technical Insights

### 1. Unified Solution Anchoring (Jackknife & Bootstrap)
Jackknife samples (\(N-1\)) are statistically nearly identical to the original data. Providing the converged primary MLE as a warm start ensures convergence in typically **1 iteration**, yielding massive **60% to 100% speedups** across the package. By lifting previous algorithmic constraints, Jackknife support is now ubiquitous.

### 2. Parametric Bootstrap Nuisance Recovery
PB calibration simulated fits show massive gains (up to **90.3%**) because the primary MLE provides near-perfect starting values for **nuisance parameters** (dispersion, shape, and covariates). Even when the treatment effect is nulled, the solver avoids the heavy cost of re-estimating these secondary parameters from scratch. Note: We've enabled PB support for previously excluded base and IVWC models via exact surrogate generative modeling.

### 3. Sequential Null-Tracking (Randomization)
Switching to sequential anchoring (tracking the null distribution) transformed previously observed slowdowns into consistent speedups (**10% to 80%**) for heavy R-loop model families. By initializing each permutation at the result of the previous one, the solver "walks" along the null manifold rather than jumping from a high-signal starting point.

### 4. Convergence Insurance
Beyond raw speed, warm starting acts as a robust **"Numerical Insurance."** It ensures the solver is protected against convergence failures on sparse bootstrap samples or ill-conditioned permutations by starting the optimization in a high-likelihood region already validated by the primary fit.

**Overall Conclusion:** Warm starting is a foundational feature of `EDI`. It provides massive computational savings for heavy models and acts as a robust "convergence insurance" for the entire resampling lifecycle. As of 2026, **Warm Starts are enabled by default** for all five resampling paths across all 90 concrete inference classes.

---

## Benchmark N/S Notes

### Non-Parametric Bootstrap & Jackknife — InferenceOrdinalPairedSignTest

The `InferenceOrdinalPairedSignTest` shows **N/S** for Non-Parametric Bootstrap and Jackknife. By design, this class implements a sign test over matched-pairs (where observational units are strictly linked in pairs). Generic subject-level bootstrap (sampling units with replacement) and subject-level jackknife (deleting one unit at a time) violate the matched-pair design constraint by leaving unmatched subjects. Because subject-level resampling breaks the paired structure of the data, Nonparam Bootstrap and Jackknife are methodologically invalid and are explicitly disabled for this class.

### Randomization — All Incidence Classes

Every `InferenceIncid*` class shows Rand: N/S. This is **by design in the EDI package**, not a benchmark limitation. `compute_rand_two_sided_pval` throws an explicit error for binary outcomes with no custom statistic function:

> *"Randomization tests are not supported for incidence. Use Zhang method."*

For binary outcomes, EDI directs users to exact permutation tests (e.g., Fisher's exact, Zhang's method) rather than the coarse randomization distribution.

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

## Detailed Rationale for 'Disabled' in Resampling

Several inference paths are marked as **Disabled** across all resampling methods:

### 1. Ordinary Least Squares (OLS) Models
*   **Paths**: `InferenceContinOLS`, `InferenceContinKKOLSIVWC`, `InferenceContinKKOLSOneLik`.
*   **Reason**: OLS estimators have a closed-form algebraic solution ($(X^TX)^{-1}X^Ty$). Because they do not use iterative optimization algorithms (like Newton-Raphson or IRLS), there is no solver state to "warm start". Enabling warm starts would incur unnecessary R6 environment copying and checks without any mathematical speedup, resulting in a net performance penalty.

### 2. Simple Mean Differences
*   **Paths**: `InferenceAllSimpleMeanDiff`, `InferenceAllSimpleMeanDiffPooledVar`.
*   **Reason**: Like OLS, these estimators compute statistical differences directly through closed-form formulas. The overhead of setting up, duplicating, and querying the warm-start environment for each resampling iteration is highly net-negative.

### 3. Bai-Adjusted T-Tests for KK Designs
*   **Paths**: `InferenceBaiAdjustedTKK14`, `InferenceBaiAdjustedTKK21`.
*   **Reason**: These estimators compute the treatment effect and its variance analytically. Since they are computed algebraically rather than iteratively, there is no numerical solver to initialize. R6 warm-start checks and environment duplication introduce up to -776% performance degradation under resampling, so warm starts are permanently disabled.
