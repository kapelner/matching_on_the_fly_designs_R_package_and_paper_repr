# Closed-Form Firth-Gradient Audit For Likelihood Paths

## Scope And Decision Rule

This audit covers the package’s likelihood-backed inference paths, grouped by the
**actual likelihood engine** they use rather than by every wrapper class that
calls the same score/Hessian code.

Operationally, I audited the paths that participate in the package’s
likelihood-backed testing/inversion framework, i.e. the classes/files that
implement `get_likelihood_test_spec()`. That is the relevant surface for the
question “can LBFGS continue uninterrupted under a Firth penalty gradient?”.

This means the audit is complete for the current **likelihood-test paths**. It
does **not** claim to cover every asymptotic model in the package, because some
classes use Wald/GEE/pass-through logic without a likelihood-test spec.

The expanded table below lists **every concrete class** on that
likelihood-backed surface. It excludes:

- abstract base classes
- alias names that point to the same concrete class object
- non-likelihood IVWC/GEE/Wald/pass-through classes without `get_likelihood_test_spec()`

The question audited is:

> Does this likelihood path admit a **closed-form Firth / Jeffreys penalty
> gradient** that is realistic enough that the current **L-BFGS** workflow could
> continue without switching to a derivative-free optimizer?

I use three labels:

- **Yes**: a closed-form gradient is standard or straightforward enough that an
  analytic adjusted score / Jeffreys-penalty gradient is realistic.
- **Borderline**: a closed-form gradient exists in principle, but would require
  a bespoke derivation that is algebraically heavy enough that I would not count
  it as a clean “LBFGS continues uninterrupted” path without dedicated work.
- **No**: the path is mixture-, latent-, quadrature-, copula-, or custom
  combined-likelihood enough that a practical closed-form Firth gradient is not
  a credible generic implementation target.

This is a **practical engineering audit**, not a statement about abstract
mathematical existence for arbitrary symbolic differentiation.

## Audit Table

### Incidence

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceIncidLogRegr` | `fast_logistic_regression_cpp` | **Yes** | Standard Bernoulli-logit GLM; classical Firth setting with realistic analytic adjusted score. |
| `InferenceIncidProbitRegr` | `fast_ordinal_probit_regression_cpp` in the 2-category case | **Yes** | Smooth fixed-effects binary probit likelihood; closed-form Jeffreys/Firth gradient is realistic in the same sense as other simple binary GLMs. |
| `InferenceIncidModifiedPoisson` | `fast_poisson_regression_cpp` | **Yes** | Same Poisson likelihood engine as ordinary Poisson regression. |
| `InferenceIncidKKModifiedPoisson` | `fast_poisson_regression_cpp` | **Yes** | Same Poisson companion likelihood engine as above. |
| `InferenceIncidLogBinomial` | `fast_log_binomial_regression_cpp` | **Yes** | Smooth binomial likelihood with noncanonical link; analytic Firth gradient is still realistic. |
| `InferenceIncidBinomialIdentityRiskDiff` | `fast_identity_binomial_regression_cpp` | **Yes** | Smooth binomial likelihood with fixed dispersion; less pleasant algebra than logit, but still tractable. |
| `InferenceIncidKKClogitOneLik` | stacked conditional-logistic + reservoir-logistic path via `fast_logistic_regression_with_var_cpp` | **Borderline** | Smooth logistic-based combined likelihood, but not the ordinary unconditional logit Firth case. |
| `InferenceIncidKKGLMM` | `fast_logistic_glmm_cpp` | **No** | GH-quadrature integrated random-effects likelihood with variance parameter. |
| `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm_cpp` | **No** | Hybrid of conditional logit and quadrature GLMM pieces with shared coefficients. |

### Count

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceCountPoisson` | `fast_poisson_regression_cpp` | **Yes** | Canonical Poisson GLM with explicit score/Hessian. |
| `InferenceCountRobustPoisson` | `fast_poisson_regression_cpp` | **Yes** | Reported estimator is robust/sandwich, but the likelihood-test path is still plain Poisson. |
| `InferenceCountQuasiPoisson` | `fast_poisson_regression_cpp` | **Yes** | Reported estimator is quasi, but the likelihood-test path is still plain Poisson. |
| `InferenceCountNegBin` | `fast_neg_bin_cpp`, `fast_neg_bin_with_var_cpp` | **Borderline** | Smooth full likelihood, but dispersion-parameter derivatives make the Firth gradient bespoke and polygamma-heavy. |
| `InferenceCountZeroInflatedPoisson` | `fast_zero_augmented_poisson_cpp` | **No** | Mixture likelihood with count and inflation blocks; practical analytic Firth gradient is not a clean target. |
| `InferenceCountHurdlePoisson` | `fast_zero_augmented_poisson_cpp` | **No** | Truncation/mixture structure makes the Jeffreys penalty gradient highly bespoke. |
| `InferenceCountZeroInflatedNegBin` | `fast_zinb_cpp` | **No** | Mixture plus dispersion parameter makes closed-form Firth support impractical. |
| `InferenceCountHurdleNegBin` | `fast_hurdle_negbin_cpp` | **No** | Same issue as above with additional hurdle/truncation structure. |
| `InferenceCountKKGLMM` | `fast_poisson_glmm_cpp` | **No** | Random-effects quadrature likelihood; not a practical closed-form Firth path. |
| `InferenceCountKKHurdlePoissonOneLik` | `fast_hurdle_poisson_glmm_cpp` | **No** | Truncated count + random effects + quadrature is too structurally complex. |
| `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined_with_var_cpp` | **No** | Hybrid conditional-plus-marginal likelihood; combined information penalty is bespoke. |

### Continuous

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceContinKKOLSOneLik` | `fast_ols_with_var_cpp` | **Yes** | Gaussian likelihood has explicit information and simple matrix derivatives. |
| `InferenceContinKKGLMM` | `fast_gaussian_lmm_cpp` | **Borderline** | Gaussian structure helps, but variance-component derivatives mean this is no longer a simple drop-in Firth path. |

### Ordinal

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceOrdinalPropOddsRegr` | `fast_ordinal_regression_cpp` | **Borderline** | Fixed-effects ordinal likelihood with thresholds; analytic path exists but needs dedicated derivation. |
| `InferenceOrdinalOrderedProbitRegr` | `fast_ordinal_probit_regression_cpp` | **Borderline** | Same threshold issue as above, with link-specific derivation for probit. |
| `InferenceOrdinalCauchitRegr` | `fast_ordinal_cauchit_regression_cpp` | **Borderline** | Smooth ordinal likelihood, but link-specific adjusted-score derivation is needed. |
| `InferenceOrdinalCloglogRegr` | `fast_ordinal_cloglog_regression_cpp` | **Borderline** | Same as above for the cloglog link. |
| `InferenceOrdinalKKGLMM` | `fast_ordinal_glmm_cpp` | **No** | Cutpoints plus quadrature plus variance parameter make this too bespoke. |

### Proportion

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferencePropBetaRegr` | `fast_beta_regression_cpp` | **Borderline** | Smooth likelihood, but mean-plus-precision structure makes the Jeffreys adjustment nontrivial. |
| `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta_cpp` | **No** | Three-component mixture; not a realistic clean analytic Firth target. |
| `InferencePropKKGLMM` | `fast_logistic_glmm_cpp` | **No** | Same engine and same quadrature/variance-parameter issue. |

### Survival

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceSurvivalCoxPHRegr` | `fast_coxph_regression_cpp` | **Borderline** | Cox bias reduction is plausible, but risk-set derivatives make this a bespoke implementation rather than a plug-in GLM case. |
| `InferenceSurvivalStratCoxPHRegr` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same Cox issue, with extra stratification structure. |
| `InferenceSurvivalKKLWACoxOneLik` | `fast_coxph_regression_cpp` | **Borderline** | Uses Cox partial likelihood over combined data; analytic bias reduction is plausible but still bespoke. |
| `InferenceSurvivalKKStratCoxOneLik` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same as above with strata-specific risk sets. |
| `InferenceSurvivalWeibullRegr` | `fast_weibull_regression_cpp` | **Borderline** | Smooth parametric survival likelihood, but materially more bespoke than the GLM cases. |
| `InferenceSurvivalDepCensTransformRegr` | `fast_dep_cens_transform_optim_cpp` | **No** | Bespoke transformation likelihood with coupled event/censoring parameter blocks. |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty_cpp` | **No** | Frailty integration and variance parameter make analytic Firth support impractical. |
| `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_clayton_weibull_aft_optim_cpp` | **No** | Copula dependence parameter plus Weibull margins plus combined design structure. |

## Concrete Conclusions

### Clean “yes” paths

These are the paths where I would expect a closed-form Firth gradient to be a
realistic extension that preserves the L-BFGS workflow:

- Bernoulli logit GLM
- Poisson GLM
- log-binomial GLM
- binomial identity-link GLM
- Gaussian linear model

These are the best first targets if the goal is to keep the present optimizer
stack and avoid numerical differentiation of the Jeffreys penalty.

### Borderline paths

These paths are analytically smooth enough that a closed-form Firth gradient is
not impossible, but I would not treat them as “LBFGS just keeps going” without
substantial model-specific derivation work:

- negative-binomial regression
- beta regression
- fixed-effects ordinal cumulative-link models, including ordered probit
- Cox / stratified Cox / LWA Cox
- Weibull regression
- Gaussian LMM
- conditional logistic matched-pair combined likelihood

In practice I would split these into two subgroups:

1. **likely worth it**:
   logit/Poisson-adjacent models, some ordinal models, maybe Cox
2. **probably not worth it early**:
   beta, negative binomial with dispersion, Gaussian LMM

### Clear “no” paths

These paths are too mixture- or latent-structure-heavy for a practical generic
closed-form Firth gradient:

- zero-inflated / hurdle count models
- zero/one-inflated beta
- all quadrature GLMM paths
- hurdle Poisson GLMM combined likelihood
- cPoisson combined likelihood
- clogit + GLMM hybrid
- frailty models
- copula models
- dependent-censoring transformation model

For these, a Firth implementation would either become:

- a bespoke research project per family, or
- a numerical penalty-gradient approximation, which defeats the goal of letting
  L-BFGS continue cleanly.

## Recommendation

If the package wants Firth support while preserving the current L-BFGS-based
architecture, I would limit the first implementation set to:

1. `InferenceIncidLogRegr`
2. `InferenceCountPoisson`
3. `InferenceIncidModifiedPoisson` and the Poisson companion-likelihood paths
4. `InferenceIncidLogBinomial`
5. `InferenceIncidBinomialIdentityRiskDiff`
6. `InferenceContinKKOLSOneLik`

After that, the next tier worth evaluating would be:

7. fixed-effects ordinal cumulative-link models, including `InferenceOrdinalOrderedProbitRegr`
8. Cox / stratified Cox
9. `InferenceIncidKKClogitOneLik`
10. maybe negative binomial

I would not plan on package-wide Firth support across the quadrature, mixture,
copula, frailty, and custom combined-likelihood engines if the requirement is
“closed-form gradient so L-BFGS continues uninterrupted.”
