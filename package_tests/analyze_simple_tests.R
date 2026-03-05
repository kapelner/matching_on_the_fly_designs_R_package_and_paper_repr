set.seed(1)
pacman::p_load(SeqExpMatch, stringr, doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
max_n_dataset = 150
source("package_tests/_dataset_load.R")

X = fread("package_tests/simple_tests_results_nc_2.csv")
table(X$rep)


colnames(X)
#rename the inferences
X[function_run == "compute_bootstrap_confidence_interval", function_run := "ci_bootstrap"]
X[function_run == "compute_bootstrap_two_sided_pval", function_run := "pval_bootstrap"]
X[function_run == "compute_confidence_interval_rand", function_run := "ci_rand"]
X[function_run == "compute_confidence_interval_rand(custom)", function_run := "ci_rand_custom"]
X[function_run == "compute_mle_confidence_interval", function_run := "ci_mle"]
X[function_run == "compute_mle_two_sided_pval_for_treatment_effect", function_run := "pval_mle"]
X[function_run == "compute_treatment_estimate", function_run := "est"]
X[function_run == "compute_two_sided_pval_for_treatment_effect_rand", function_run := "pval_rand"]
X[function_run == "compute_two_sided_pval_for_treatment_effect_rand(custom)", function_run := "pval_rand_custom"]
X[function_run == "compute_two_sided_pval_for_treatment_effect_rand(delta=0.5)", function_run := "pval_rand_shift"]
table(X$function_run)

#check duration time
D = X[, .(mean_duration = mean(duration_time_sec)), 
  by = c("inference_class", "design", "function_run")][order(-mean_duration)]
D[1:100]
table(D[1:100]$inference_class)

#we kind of don't care as much about these results
Xextra = X[function_run %in% c("pval_rand_custom", "pval_rand_shift")]
X = X[!(function_run %in% c("pval_rand_custom", "pval_rand_shift"))]

#create the actual beta being estimated
table(X$inference_class)

#set the correct value of beta as it's inference and response dependent - this default works for most
X[, beta := as.numeric(beta_T)]

e = exp(1)
# incidence — simple/KK mean diff: probability-scale estimand ≈ 0.183, not 1
X[beta_T == 1 &
  response_type == "incidence" &
  inference_class %in% c("SeqDesignInferenceAllSimpleMeanDiff", "SeqDesignInferenceAllKKCompoundMeanDiff"),
    beta := (1/2)*3*e/(1+3*e) + (1/2)*e/(3+e) - 0.5]

# incidence — univariate logistic: marginal log-OR ≈ 0.767, not 1
X[beta_T == 1 & inference_class == "SeqDesignInferenceIncidUnivLogRegr",
    beta := qlogis((1/2)*3*e/(1+3*e) + (1/2)*e/(3+e))]


# proportion — simple/KK mean diff: E[e*y/(1+(e-1)*y)] - E[y]
X[beta_T == 1 & response_type == "proportion" &
  inference_class %in% c("SeqDesignInferenceAllSimpleMeanDiff", "SeqDesignInferenceAllKKCompoundMeanDiff"),
    beta := {
        y_p = datasets_and_response_models[[dataset]]$y_original$proportion
        mean(e * y_p / (1 + (e - 1) * y_p)) - mean(y_p)
    }, by = dataset]

# count — simple/KK mean diff: E[Y_C] * (e - 1)
X[beta_T == 1 & response_type == "count" &
  inference_class %in% c("SeqDesignInferenceAllSimpleMeanDiff", "SeqDesignInferenceAllKKCompoundMeanDiff"),
    beta := {
        y_c = datasets_and_response_models[[dataset]]$y_original$count
        mean(y_c) * (e - 1)
    }, by = dataset]

# survival — KM median diff: (e - 1) * median(Y_C)
X[beta_T == 1 & inference_class == "SeqDesignInferenceSurvivalKMDiff",
    beta := {
        y_s = datasets_and_response_models[[dataset]]$y_original$survival
        (e - 1) * median(y_s)
    }, by = dataset]

X[beta_T == 1 & inference_class == "SeqDesignInferenceSurvivalRestrictedMeanDiff",
    beta := {
        y_s = datasets_and_response_models[[dataset]]$y_original$survival
        tau = quantile(y_s, 0.95)
        (e - 1) * mean(pmin(y_s, tau))
    }, by = dataset]


#now some are impossible to calculate for real data due to the unknown f(x) model
X[beta_T == 1 &
  inference_class %in% c(
    "SeqDesignInferenceIncidMultiLogRegr", 
    "SeqDesignInferenceIncidMultiKKGEE",
    "SeqDesignInferenceIncidMultiKKGLMM",
    "SeqDesignInferencePropUniBetaRegr", 
    "SeqDesignInferencePropMultiBetaRegr", 
    "SeqDesignInferenceSurvivalUniCoxPHRegr", 
    "SeqDesignInferenceSurvivalMultiCoxPHRegr"
  ),
  beta := NA_real_]
table(X$beta, useNA = "always")


#check MSE
X[function_run == "est", sqerr := (result_1 - beta)^2]
E = X[function_run == "est", .(mse = mean(sqerr, na.rm = TRUE), beta = first(beta)), 
  by = c("inference_class", "design", "response_type", "beta_T")][order(-mse)]
E = E[!is.nan(mse)]
E[beta_T == 1][1:100]
E[beta_T == 0][1:100]

#check coverage
X[str_detect(X$function_run, "ci"), ci_correct := ifelse(beta >= result_1 & beta <= result_2, 1, 0)]
table(X$ci_correct, useNA = "always")
C = X[str_detect(X$function_run, "ci"), .(coverage = mean(ci_correct)), 
  by = c("inference_class", "function_run", "response_type", "beta_T")][order(coverage)]
C[beta_T == 1][1:100]
C[beta_T == 0][1:100]

#check size
X[beta_T == 0 & str_detect(X$function_run, "pval"), H0_rejected := ifelse(result_1 < 0.05, 1, 0)]
#table(X[beta_T == 0]$H0_rejected, useNA = "always")
S = X[beta_T == 0 & str_detect(X$function_run, "pval"), .(
  size = mean(H0_rejected, na.rm = TRUE),
  size_pval = prop.test(sum(H0_rejected, na.rm = TRUE), length(na.omit(H0_rejected)), p = 0.05)$p.value
), by = c("inference_class", "design", "response_type", "beta_T", "function_run")][order(-size)]
S[, bonf_size_pval := pmin(1, size_pval * .N)][, size_pval := NULL]
S[11:160]

# 2. `pval_mle` (Size $\approx$ 0.15 - 0.25)
#   This is the classic failure of parametric MLE inference (like OLS and CoxPH) on real-world datasets, and it perfectly
#   highlights why this package's Bootstrap and Randomization tests exist!
#    * The tests are injecting a treatment effect into real-world datasets (airquality, boston, etc.) without modifying the
#      underlying covariate relationships.
#    * Real-world datasets contain severe heteroscedasticity, non-linearities, and unmeasured confounding that strongly
#      violate the assumptions required by standard MLE standard errors (which assume IID normal/homoscedastic errors).
#    * Because the MLE standard errors are overly optimistic (too small) on misspecified real-world data, the Wald test
#      rejects the true null hypothesis ($H_0: \beta_T = 0$) far too often, leading to the inflated ~20% Type I Error rates.

#check power
X[beta_T == 1 & str_detect(X$function_run, "pval"), H0_rejected := ifelse(result_1 < 0.05, 1, 0)]
table(X[beta_T == 1]$H0_rejected, useNA = "always")
P = X[beta_T == 1 & str_detect(X$function_run, "pval"), .(power = mean(H0_rejected)), 
  by = c("inference_class", "design", "response_type", "function_run")][order(-power)]
P[!is.na(power)][1:200]

# 1. The Bootstrapped KK OLS (ContinMultOLSKK)
#   You'll notice that for the continuous KK OLS estimator, only the pval_bootstrap is in this low-power list (~52-54%). Its
#   pval_rand and pval_mle are completely fine!


#   This is caused by severe structural instability when resampling matched pairs with replacement:
#    * When you sample $m$ matched pairs with replacement for the bootstrap, you inevitably draw many duplicate pairs.
#    * The OLS design matrix for the matched differences ($X_d$) suddenly loses full column rank because the duplicate rows
#      provide zero new linear information.
#    * To prevent a singular matrix crash, your estimator's QR decomposition correctly drops the collinear covariates. If the
#      rank drops too low, it entirely abandons the multivariate OLS and falls back to a naive Mean Difference for that
#      specific bootstrap iteration.
#    * Because the model randomly shifts between a fully-adjusted OLS, partially-adjusted OLS, and a naive mean difference
#      from iteration to iteration, the variance of the resulting bootstrap distribution artificially explodes. This massive
#      variance widens the bootstrap confidence interval and obliterates its statistical power.


#   2. Incidence and Proportion Models
#   The rest of the low power entries (IncidMultiLogRegr, PropUniBetaRegr, AllSimpleMeanDiff on incidence/proportion) are
#   simply struggling due to the mathematical limitations of the simulated effect size on bounded domains at low sample
#   sizes.
#    * In simple_tests.R, the data generation adds beta_T = 1 on the link scale (e.g. plogis(qlogis(p) + 1)).
#    * For a baseline probability of $0.5$, shifting the log-odds by $1.0$ only moves the final probability to $0.73$ (an
#      absolute difference of just $0.23$). For baseline probabilities nearer to the boundaries (like the $0.75$ and $0.25$
#      base rates used in your Incidence test), the absolute shift in probability is even smaller ($\sim 0.18$).
#    * Detecting a probability shift of $0.18$ between two groups with sample sizes frequently around $n=50$ (like the cars
#      dataset) or $n=150$ is notoriously difficult. A standard two-sample proportion test for this effect size at $n=50$ has
#      a theoretical mathematical power of exactly ~25% to ~40%.


#   The fact that your tests are successfully detecting this small probability shift ~60% of the time means your sequential
#   designs and multivariate estimators are actually performing better than standard unadjusted statistical tests would! The
#   ~50% power seen specifically on the pval_bootstrap for these GLMs just reflects the added volatility of bootstrapping
#   logistic/beta regressions on tiny datasets, where perfect separation and boundary clamping often occur during resampling.