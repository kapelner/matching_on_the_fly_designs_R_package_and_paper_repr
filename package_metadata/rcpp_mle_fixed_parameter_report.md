# Fixed-Parameter Optimization for Rcpp MLE Implementations

This report covers the package's internal Rcpp MLE or M-estimation fitters where
we optimize parameters ourselves. The requested constraint is an equality
constraint: some entries of the parameter vector are fixed before optimization,
and only the remaining entries are estimated.

## Recommended General Design

Use a reduced-parameter adapter, not a new constrained optimizer.

For every fitter with full parameter vector `theta`, accept optional arguments:

- `fixed_idx`: one-based integer indices of parameters fixed by the caller.
- `fixed_values`: numeric values for those entries.
- `start_params`: full-length or free-only start values; full-length is safer.

Internally:

1. Build `free_idx = setdiff(seq_along(theta), fixed_idx)`.
2. Optimize only `theta_free`.
3. On each objective/score evaluation, reconstruct the full vector:
   `theta[fixed_idx] = fixed_values`, `theta[free_idx] = theta_free`.
4. Return the full parameter vector, plus `fixed_idx`, `free_idx`, and a
   convergence flag.
5. For standard errors, compute the observed information or Hessian for the full
   vector at the constrained MLE, then invert only the free-free block. Fixed
   parameters should have variance 0 if treated as known constants, or `NA` if
   downstream code should not report an SE for fixed entries. Cross-covariances
   involving fixed parameters should be 0 or `NA` consistently.

This keeps the likelihood code unchanged and avoids introducing box-constraint
machinery for a problem that is just parameter masking.

## Cross-Cutting Difficulty

The work is mostly mechanical for likelihood functors that already expose
`operator()` and `hessian()`: wrap the functor with a free-to-full parameter map
and pass the reduced gradient to LBFGS/Newton. The harder part is not the
optimizer; it is making every return shape and every downstream standard-error
consumer understand fixed parameters.

The highest-risk cases are models with internal parameter constraints or
identified-by-construction parameterizations, especially ordinal thresholds,
stereotype-logit category scores, and variance components.

## Implementation Report

| Rcpp implementation | Used by inference paths | Current optimizer | Difficulty | Why / required changes |
|---|---|---:|---:|---|
| `fast_logistic_regression_cpp`, `fast_logistic_regression_weighted_cpp`, `fast_logistic_regression_with_var_cpp` | `LogRegr`, GComp incidence/proportion, `KKClogit*`, CPoisson matched pieces | IRLS | Easy | Treat fixed coefficients as an offset: `eta_fixed = X_fixed beta_fixed`, run IRLS on `X_free` with adjusted linear predictor. Variance is inverse free-free Fisher block. |
| `fast_poisson_regression_cpp`, weighted, quasi, robust-Poisson working paths | `Poisson`, `QuasiPoisson`, `ModifiedPoisson`, `KKModifiedPoisson`, `RobustPoisson`, CPoisson reservoir fallback | IRLS | Easy | Same offset pattern as logistic. For quasi-Poisson, dispersion is computed after the constrained fit and the covariance is the free-free inverse information scaled by dispersion. |
| `fast_ols_cpp`, `fast_ols_with_var_cpp` | `OLS`, `RiskDiff`, linear working estimators | closed-form QR solve | Easy | Move fixed terms to the response: `y_adj = y - X_fixed beta_fixed`, then QR-solve on `X_free`. This is the simplest constrained path. |
| `fast_robust_regression_cpp` | `RobustRegr`, `KKRobustRegr*` | robust IRLS | Easy | Same `y_adj`/`X_free` approach as OLS inside each weighted least-squares step. Need to compute residuals using full `X theta`. |
| `fast_log_binomial_regression_cpp`, `fast_identity_binomial_regression_cpp` | `LogBinomial`, `BinomialIdentityRiskDiff` | constrained Fisher scoring with step-halving | Moderate | Use a fixed linear-predictor offset and update only free coefficients. The feasibility checks for `mu` must include the fixed offset. Step-halving remains valid but should evaluate full constrained likelihood. |
| `fast_beta_regression_cpp`, `fast_beta_regression_with_var_cpp` | `BetaRegr` | LBFGS | Easy | Parameter vector is `[beta, log_phi]`. Add a free-parameter LBFGS wrapper around the existing beta likelihood. Return full coefficients and invert the free-free Hessian. Fixing `log_phi` is straightforward. |
| `fast_zero_one_inflated_beta_cpp` | `ZeroOneInflatedBetaRegr` | LBFGS | Easy to Moderate | The likelihood already has a single full vector for mean, zero inflation, one inflation, and beta precision pieces. The generic LBFGS masking adapter works. Risk is documenting the parameter index layout clearly enough for callers. |
| `fast_zero_augmented_poisson_cpp` | `ZeroInflatedPoisson`, `HurdlePoisson` | LBFGS | Easy to Moderate | Parameter vector is split into conditional and zero/hurdle components. Masking works directly. Need to preserve support for `estimate_only=TRUE` and ensure fixed parameters are respected by both conditional and zero/hurdle design matrices. |
| `fast_neg_bin_cpp`, `fast_neg_bin_with_var_cpp` | `NegBin`, CPoisson IVWC reservoir | LBFGS | Easy | Parameter vector is `[beta, log_theta]`. Generic LBFGS wrapper is enough. Variance uses free-free Hessian. Fixing overdispersion is straightforward. |
| `fast_hurdle_negbin_cpp`, `fast_hurdle_negbin_with_var_cpp` | `HurdleNegBin` | logistic IRLS plus LBFGS truncated-NB count component | Moderate | There are two model components. The hurdle/logistic part can use the logistic constrained path; the positive-count truncated NB part can use LBFGS masking. Need an API that can fix parameters separately for hurdle and count components, or a documented combined index layout. |
| `fast_weibull_regression_cpp` | `WeibullRegr` | LBFGS | Easy | Parameter vector includes AFT regression parameters and log-scale/log-shape component. Generic LBFGS masking applies. Need clear mapping for the non-beta survival parameter. |
| `fast_coxph_regression_cpp` | `CoxPHRegr` | Newton-Raphson | Moderate | The Cox partial likelihood has no intercept and updates `beta` via Newton steps. Fixed coefficients can be handled as a risk-score offset with Newton on `X_free`. Need careful updates to risk-set sums so `exp(X_fixed beta_fixed + X_free beta_free)` is used everywhere. |
| `fast_dep_cens_transform_optim_cpp` | `DepCensTransformRegr` | LBFGS | Easy | Existing likelihood functor has full `operator()` and `hessian()`. Generic LBFGS masking is enough. The recent vcov-hardening logic should operate on the free-free Hessian. |
| `fast_clayton_weibull_aft_optim_cpp` | `KKClaytonCopulaIVWC`, `KKClaytonCopulaOneLik` | LBFGS | Easy to Moderate | Generic masking works for the full vector, but this model has dependence and Weibull nuisance parameters. Need a stable documented index map and tests for fixing the dependence parameter versus fixing treatment/covariate coefficients. |
| `fast_clogit_plus_glmm_cpp` | `KKClogitPlusGLMMIVWC`, `KKClogitPlusGLMMOneLik` | LBFGS | Easy to Moderate | The combined likelihood already evaluates a full parameter vector. Masking works. The main risk is treatment-column indexing across discordant and concordant components and preserving the `estimate_only=TRUE` hessian-free path. |
| `fast_cpoisson_combined_with_var_cpp` | `KKCPoissonOneLik` | Newton-Raphson using analytic information | Moderate | The score/information are already explicit. For fixed parameters, solve only `info_free_free * step = score_free`, update free entries, and keep fixed entries constant. Variance is inverse free-free information. This is clean but needs custom Newton edits. |
| `fast_gaussian_lmm_cpp` | `KKGLMM` continuous fast path | LBFGS | Moderate | Generic LBFGS masking works, but variance components need constraints and interpretation. Fixing fixed-effect betas is easy; fixing log-variance components is also possible but should be separately tested because boundary behavior and Hessian conditioning are more fragile. |
| `fast_ordinal_regression_cpp` | `PropOddsRegr`, ordinal GComp | Newton-Raphson with finite-difference derivatives and line search | Moderate to Hard | Parameter vector is `[thresholds, beta]` with ordered thresholds. Masking works mechanically, but fixed thresholds can violate ordering or make the free threshold constraints awkward. We need validation that fixed thresholds are strictly ordered relative to neighboring free thresholds. |
| `fast_ordinal_probit_regression_cpp` | `OrderedProbitRegr` | Newton-Raphson with finite-difference derivatives and line search | Moderate to Hard | Same threshold-ordering issue as proportional odds, with probit link. Generic free-vector Newton is possible, but constraint validation is required. |
| `fast_ordinal_cloglog_regression_cpp` | `CloglogRegr` | Newton-Raphson with finite-difference derivatives and line search | Moderate to Hard | Same threshold-ordering issue as proportional odds. |
| `fast_ordinal_cauchit_regression_cpp` | `CauchitRegr` | Newton-Raphson with finite-difference derivatives and line search | Moderate to Hard | Same threshold-ordering issue as proportional odds; the cauchit tail behavior may make bad fixed threshold values especially unstable. |
| `fast_adjacent_category_logit_cpp` | Rcpp implementation present; not currently exported as an inference path in metadata | Newton-like finite-difference / full-vector updates | Moderate | The full vector has category intercept/slope structure. Masking is feasible but needs a documented parameter layout. Variance is free-free Hessian. |
| `fast_continuation_ratio_regression_cpp` | Rcpp implementation present; not currently exported as an inference path in metadata | IRLS/Newton-style iterative fit | Moderate | Often representable as augmented binary regressions. Fixed parameters can be handled either through offsets in each augmented logistic block or by full-vector masking. Need care because one logical coefficient may appear in repeated augmented rows. |
| `fast_stereotype_logit_cpp` | Rcpp implementation present; not currently exported as an inference path in metadata | nested/profile optimization | Hard | This is the hardest. It has nuisance/category score parameters with identification constraints and a profiling helper already exists for fixed `beta`. General fixed arbitrary parameters would need custom handling for the score-parameter transformation and nested optimizer. |

## API Shape

Recommended exported signatures should add optional arguments at the end to avoid
breaking existing R calls:

```r
fast_model_cpp(..., fixed_idx = integer(0), fixed_values = numeric(0))
fast_model_with_var_cpp(..., fixed_idx = integer(0), fixed_values = numeric(0))
```

For multi-component models, use named lists on the R side if a single global
index layout is too error-prone:

```r
fixed = list(
  cond = c("(Intercept)" = 0),
  zi = c("w" = 0),
  dispersion = c(log_theta = log(5))
)
```

The R wrapper can translate names to global C++ indices before calling Rcpp.
That keeps the C++ interface fast and simple while making the public R API safer.

## Variance and Inference Rules

For a constrained MLE where fixed parameters are treated as known constants:

- Optimize and report only free-parameter uncertainty.
- Compute `I_ff`, the observed information submatrix for free parameters.
- Use `solve(I_ff)` for the covariance of free estimates.
- Return a full covariance matrix with fixed rows/columns set to 0 or `NA`.
  Prefer `NA` if callers might accidentally interpret a fixed parameter as
  estimated with zero uncertainty.
- If the treatment parameter is fixed, asymptotic CI/p-value methods for that
  parameter should return `NA` or a clear non-estimable status rather than a zero
  SE Wald result.

## Suggested Rollout

1. Add shared C++ helpers:
   - validate fixed indices and values;
   - build `free_idx`;
   - expand free parameters to full parameters;
   - subset full gradients/Hessians to free entries;
   - expand free covariance back to full shape.
2. Implement first for the generic LBFGS functor models:
   `beta`, `negative binomial`, `Weibull`, survival optimizers,
   zero-augmented Poisson, zero-one-inflated beta, and PlusGLMM.
3. Implement IRLS offset models:
   logistic, Poisson, OLS, robust, log-binomial, identity-binomial.
4. Implement custom Newton models:
   Cox, CPoisson combined, ordinal cumulative-link models.
5. Leave stereotype-logit for last because its identification constraints make
   arbitrary fixed-parameter support substantially more delicate.

## Overall Estimate

This is feasible without replacing optimizers.

- Easy models: about half a day to one day for a shared helper plus the main
  LBFGS/IRLS families.
- Moderate models: one to three days for Cox, CPoisson combined, hurdle NB,
  Gaussian LMM, and ordinal cumulative links with tests.
- Hard model: stereotype logit could take another one to two days depending on
  how general the fixed-parameter API must be.

The most practical path is to support fixed treatment/covariate coefficients
first, then decide whether fixing thresholds, dispersion, shape, and dependence
parameters is required. Arbitrary fixed nuisance parameters are possible, but
they increase testing and documentation burden more than they increase optimizer
complexity.
