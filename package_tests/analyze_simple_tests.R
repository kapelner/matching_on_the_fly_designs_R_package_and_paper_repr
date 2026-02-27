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

#create the actual beta being estimated
table(X$inference_class)
#good default


D = X[, .(mean_duration = mean(duration_time_sec)), 
  by = c("inference_class", "design", "function_run")][order(-mean_duration)]
table(D[1:100]$inference_class)
table(X$function_run)


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
  inference_class %in% c("SeqDesignInferenceIncidMultiLogRegr", "SeqDesignInferencePropUniBetaRegr", "SeqDesignInferencePropMultiBetaRegr", "SeqDesignInferenceSurvivalUniCoxPHRegr", "SeqDesignInferenceSurvivalMultiCoxPHRegr"),
    beta := NA_real_]
table(X$beta, useNA = "always")


#check MSE
X[function_run == "est", sqerr := (result_1 - beta)^2]
E = X[function_run == "est", .(mse = mean(sqerr), beta = first(beta)), 
  by = c("inference_class", "design", "response_type", "beta_T")][order(-mse)]
E[beta_T == 1][1:100]
E[beta_T == 0][1:100]

#check coverage
X[str_detect(X$function_run, "ci"), ci_correct := ifelse(beta >= result_1 & beta <= result_2, 1, 0)]
table(X$ci_correct, useNA = "always")
C = X[str_detect(X$function_run, "ci"), .(coverage = mean(ci_correct)), 
  by = c("inference_class", "function_run", "response_type", "beta_T")][order(coverage)]
C[beta_T == 1][1:100]
C[beta_T == 0][1:100]

#check power
X[str_detect(X$function_run, "pval"), pval_correct := ifelse(beta_T == 1 & result_1 < 0.05, 1, ifelse(beta_T == 0 & result_1 > 0.05, 1, 0))]
table(X$pval_correct, useNA = "always")
P = X[str_detect(X$function_run, "pval"), .(power = mean(pval_correct)), 
  by = c("inference_class", "design", "response_type", "beta_T")][order(-power)]
P[beta_T == 1][1:100]
P[beta_T == 0][1:100]

