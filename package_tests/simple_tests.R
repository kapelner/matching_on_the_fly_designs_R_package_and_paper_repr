#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.5/ SeqExpMatch/
rm(list = ls())
set.seed(1)
pacman::p_load(SeqExpMatch, doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
max_n_dataset = 150
source("package_tests/_dataset_load.R")
options(error = recover)
# options(warn=2)

Nrep = 10
NUM_CORES = 2
prob_censoring = 0.15
nsim_exact_test = 351
pval_epsilon = 0.007
test_compute_confidence_interval_rand = TRUE #too slow for now
beta_T_values = c(0, 1)
SD_NOISE = 0.1

results_file = paste0("package_tests/simple_tests_results_nc_", NUM_CORES, ".csv")
run_row_id = 0L
results_dt = data.table(
  duration_time_sec = numeric(),
  inference_class = character(),
  dataset = character(),
  design = character(),
  response_type = character(),
  function_run = character(),
  result_1 = character(),
  result_2 = character(),
  beta_T_in_confidence_interval = logical(),
  error_message = character(),
  rep = integer(),
  run_row_id = integer(),
  beta_T = numeric(),
  nsim_exact_test = integer(),
  pval_epsilon = numeric(),
  prob_censoring = numeric(),
  sd_noise = numeric(),
  num_cores = integer(),
  dataset_n_rows = integer(),
  dataset_n_cols = integer(),
  result = character(),
  status = character()
)
write_results_if_needed = function(force = FALSE){
  if (force || nrow(results_dt) > 0L){
    data.table::fwrite(results_dt, results_file)
  }
}

record_result = function(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, inference_class, function_run, result, status, duration_time_sec, error_message = NA_character_){
  result_vec = if (is.null(result)) {
    NA_character_
  } else if (length(result) == 0) {
    character(0)
  } else {
    as.character(result)
  }

  result_str = if (is.null(result)) {
    NA_character_
  } else if (length(result) == 0) {
    ""
  } else if (length(result) == 1) {
    as.character(result)
  } else {
    paste(as.character(result), collapse = " ")
  }
  result_1 = if (length(result_vec) >= 1) result_vec[1] else NA_character_
  result_2 = if (grepl("confidence_interval", function_run, fixed = TRUE) && length(result_vec) >= 2) result_vec[2] else NA_character_
  beta_T_in_confidence_interval = NA
  if (grepl("confidence_interval", function_run, fixed = TRUE) && length(result) >= 2 && all(is.finite(result[1:2]))){
    ci_lo = min(result[1:2])
    ci_hi = max(result[1:2])
    beta_T_in_confidence_interval = (beta_T >= ci_lo && beta_T <= ci_hi)
  }
  run_row_id <<- run_row_id + 1L
  results_dt <<- data.table::rbindlist(list(
    results_dt,
    data.table(
      rep = as.integer(rep_curr),
      run_row_id = run_row_id,
      duration_time_sec = duration_time_sec,
      beta_T = beta_T,
      nsim_exact_test = as.integer(nsim_exact_test),
      pval_epsilon = pval_epsilon,
      prob_censoring = prob_censoring,
      sd_noise = SD_NOISE,
      num_cores = as.integer(NUM_CORES),
      inference_class = inference_class,
      dataset = dataset_name,
      dataset_n_rows = as.integer(dataset_n_rows),
      dataset_n_cols = as.integer(dataset_n_cols),
      design = design_type,
      response_type = response_type,
      function_run = function_run,
      result_1 = result_1,
      result_2 = result_2,
      beta_T_in_confidence_interval = beta_T_in_confidence_interval,
      result = result_str,
      status = status,
      error_message = ifelse(is.null(error_message), NA_character_, as.character(error_message))
    )
  ), use.names = TRUE)
  write_results_if_needed(force = TRUE)
}

run_inference_checks = function(seq_des_inf, response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols){
  skip_slow = is(seq_des_inf, "SeqDesignInferencePropMultiBetaRegr") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiWeibullRegr") || is(seq_des_inf, "SeqDesignInferenceCountMultiNegBinRegr") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiCoxPHRegr")
  skip_ci = beta_T == 1 && (
    is(seq_des_inf, "SeqDesignInferenceIncidMultiLogRegr") ||
    is(seq_des_inf, "SeqDesignInferencePropUniBetaRegr") ||
    is(seq_des_inf, "SeqDesignInferencePropMultiBetaRegr") ||
    is(seq_des_inf, "SeqDesignInferenceSurvivalUniCoxPHRegr") ||
    is(seq_des_inf, "SeqDesignInferenceSurvivalMultiCoxPHRegr")
  )
  snap_small_numeric_to_zero = function(x, tol = sqrt(.Machine$double.eps)){
    if (is.null(x)){
      return(x)
    }
    if (is.list(x)){
      return(lapply(x, snap_small_numeric_to_zero, tol = tol))
    }
    if (is.atomic(x) && is.numeric(x)){
      x[is.finite(x) & abs(x) < tol] = 0
      return(x)
    }
    x
  }

  has_invalid_numeric = function(x){
    if (is.null(x)){
      return(FALSE)
    }
    if (is.list(x)){
      return(any(vapply(x, has_invalid_numeric, logical(1))))
    }
    if (is.atomic(x) && is.numeric(x)){
      return(any(!is.finite(x) | is.na(x) | is.nan(x)))
    }
    FALSE
  }

  is_zero_zero_confidence_interval = function(label, result){
    if (!grepl("confidence_interval", label, fixed = TRUE)){
      return(FALSE)
    }
    if (!(is.atomic(result) && is.numeric(result))){
      return(FALSE)
    }
    if (length(result) < 2){
      return(FALSE)
    }
    isTRUE(all(result[1:2] == 0))
  }

  safe_call = function(label, expr){
    start_elapsed = unname(proc.time()[["elapsed"]])
    tryCatch({
        result <- expr
        if (has_invalid_numeric(result)){
          stop("Invalid output detected (NA/NaN/Inf) in ", label)
        }
        if (is_zero_zero_confidence_interval(label, result)){
          stop("Degenerate confidence interval [0, 0] detected in ", label)
        }
        result = snap_small_numeric_to_zero(result)
        cat("  ", result, "\n")
        duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
        record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, class(seq_des_inf)[1], label, result, status = "ok", duration_time_sec = duration_time_sec)
        result
      }, error = function(e){
      is_non_fatal = grepl("not implemented", e$message, fixed = TRUE) ||
                     grepl("not enough discordant pairs", e$message, ignore.case = TRUE) ||
                     grepl("Degenerate confidence interval", e$message, fixed = TRUE) ||
                     grepl("inconsistent estimator units", e$message, ignore.case = TRUE) ||
                     grepl("Bootstrap confidence interval returned NA bounds", e$message, fixed = TRUE) ||
                     grepl("Weibull regression failed to converge", e$message, fixed = TRUE) ||
                     (grepl("NA/NaN/Inf", e$message, fixed = TRUE) &&
                      grepl("compute_bootstrap", label, fixed = TRUE) &&
                      (is(seq_des_inf, "SeqDesignInferenceIncidKKClogit") || is(seq_des_inf, "SeqDesignInferenceContinMultGLS")))
      
      if (is_non_fatal){
        message("Skipping ", label, " (non-fatal): ", e$message)
        duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
        record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, class(seq_des_inf)[1], label, NA_character_, status = "error", duration_time_sec = duration_time_sec, error_message = e$message)
      } else {
        stop(e$message)
      }
    })
  }

  message("    Calling compute_treatment_estimate()")
  safe_call("compute_treatment_estimate",
    seq_des_inf$compute_treatment_estimate())
  message("    Calling compute_mle_two_sided_pval_for_treatment_effect()")
  safe_call("compute_mle_two_sided_pval_for_treatment_effect",
            seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect())
  if (!skip_ci){
    message("    Calling compute_mle_confidence_interval()")
    safe_call("compute_mle_confidence_interval",
              seq_des_inf$compute_mle_confidence_interval(0.05))
  }
  if (!skip_slow && !skip_ci){
    message("    Calling compute_bootstrap_confidence_interval()")
    safe_call("compute_bootstrap_confidence_interval",
              seq_des_inf$compute_bootstrap_confidence_interval(B = nsim_exact_test, na.rm = TRUE))
  }
  if (!skip_slow){
    message("    Calling compute_bootstrap_two_sided_pval()")
    safe_call("compute_bootstrap_two_sided_pval",
              seq_des_inf$compute_bootstrap_two_sided_pval(B = nsim_exact_test, na.rm = TRUE))
  }
  if (!skip_slow && response_type %in% c("continuous", "survival", "proportion")){
    message("    Calling compute_two_sided_pval_for_treatment_effect_rand()")
    safe_call("compute_two_sided_pval_for_treatment_effect_rand",
            seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = nsim_exact_test, show_progress = FALSE))
    message("    Calling compute_two_sided_pval_for_treatment_effect_rand(delta=0.5)")
    use_log   = response_type == "survival"
    use_logit = response_type == "proportion"
    transform_for_rand = if (use_log) "log" else if (use_logit) "logit" else "none"
    delta_for_rand = 0.5  # on the transformed scale: log-ratio for survival, log-odds-ratio for proportion, additive for continuous
    safe_call("compute_two_sided_pval_for_treatment_effect_rand(delta=0.5)",
            seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(
              nsim_exact_test = nsim_exact_test,
              delta = delta_for_rand,
              transform_responses = transform_for_rand,
              show_progress = FALSE
            ))
  }

  if (!skip_slow && !skip_ci && test_compute_confidence_interval_rand & response_type %in% c("continuous",  "proportion",  "count")){ #,  "proportion", "survival"
    message("    Calling compute_confidence_interval_rand()")
    safe_call("compute_confidence_interval_rand",
              seq_des_inf$compute_confidence_interval_rand(nsim_exact_test = nsim_exact_test, pval_epsilon = pval_epsilon, show_progress = FALSE))
    #message("    Calling compute_confidence_interval_rand(alpha=0.01)")
    #safe_call("compute_confidence_interval_rand(alpha=0.01)",
    #          seq_des_inf$compute_confidence_interval_rand(nsim_exact_test = nsim_exact_test, pval_epsilon = pval_epsilon, alpha = 0.01))
  }
  message("    Calling set_custom_randomization_statistic_function()")
  seq_des_inf$set_custom_randomization_statistic_function(function(){
    yTs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 1]
    yCs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 0]
    (mean(yTs) - mean(yCs)) / sqrt(var(yTs) / length(yTs) + var(yCs) / length(yCs))
  })
  if (!skip_slow){
    message("    Calling compute_two_sided_pval_for_treatment_effect_rand(custom)")
    safe_call("compute_two_sided_pval_for_treatment_effect_rand(custom)",
              seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = nsim_exact_test, show_progress = FALSE))
  }
  if (!skip_slow && !skip_ci && test_compute_confidence_interval_rand & response_type %in% c("continuous")){ #, "proportion", "survival"
    message("    Calling compute_confidence_interval_rand(custom)")
    safe_call("compute_confidence_interval_rand(custom)",
              seq_des_inf$compute_confidence_interval_rand(nsim_exact_test = nsim_exact_test, pval_epsilon = pval_epsilon, show_progress = FALSE))
  }
  seq_des_inf$set_custom_randomization_statistic_function(NULL)
}

run_tests_for_response = function(response_type, design_type, dataset_name){
  inference_banner = function(inf_name){
    message(paste0("\n\n  == Inference: ", inf_name, " dataset = ", dataset_name, " beta_T = [", beta_T, "] num_cores = [", NUM_CORES, "] rep = [", rep_curr, "/", Nrep, "]\n"))
  }

  apply_treatment_effect_and_noise = function(y_t, w_t, response_type){
    eps = rnorm(1, 0, SD_NOISE)
    bt = ifelse(w_t == 1, beta_T, 0)    
    if (response_type == "continuous"){
      return(y_t + bt + eps)
    }
    if (response_type == "incidence"){
      p_base = ifelse(y_t > 0, 0.75, 0.25)
      p_t = plogis(qlogis(p_base) + bt + eps)
      return(as.numeric(stats::rbinom(1, size = 1, prob = p_t)))
    }
    if (response_type == "proportion"){
      y_clamped = pmax(.Machine$double.eps, pmin(1 - .Machine$double.eps, y_t))
      return(plogis(qlogis(y_clamped) + bt + eps))
    }
    if (response_type == "count"){
      lambda_t = pmax(.Machine$double.eps, y_t * exp(bt + eps))
      return(as.numeric(stats::rpois(1, lambda = lambda_t)))
    }
    if (response_type == "survival"){
      return(pmax(.Machine$double.eps, y_t * exp(bt + eps)))
    }
    stop("Unsupported response_type for treatment effect: ", response_type)
  }

  D = datasets_and_response_models[[dataset_name]]
  n = nrow(D$X)
  dead = rep(1, n)
  y = D$y_original[[response_type]]
  t_f = quantile(y, .95)

  if (response_type == "survival"){
    for (i in 1 : n){
      if (runif(1) < prob_censoring | y[i] >= t_f){
        y[i] = runif(1, 0, y[i])
        dead[i] = 0
      }
    }
  }

  seq_des_obj = switch(design_type,
    KK21 =         SeqDesignKK21$new(        response_type = response_type, n = n),
    KK21stepwise = SeqDesignKK21stepwise$new(response_type = response_type, n = n),
    KK14 =         SeqDesignKK14$new(        response_type = response_type, n = n),
    CRD =          SeqDesignCRD$new(         response_type = response_type, n = n),
    Efron =        SeqDesignEfron$new(       response_type = response_type, n = n),
    Atkinson =     SeqDesignAtkinson$new(    response_type = response_type, n = n),
    iBCRD =        SeqDesigniBCRD$new(       response_type = response_type, n = n),
    stop("Unsupported design_type: ", design_type)
  )
  for (t in 1 : n){
    w_t = seq_des_obj$add_subject_to_experiment_and_assign(D$X[t, ])
    y_t = apply_treatment_effect_and_noise(y[t], w_t, response_type)
    seq_des_obj$add_subject_response(t, y_t, dead[t])
  }

  is_kk_design = design_type %in% c("KK21", "KK21stepwise", "KK14")
  if (response_type == "continuous"){
    inference_banner("SeqDesignInferenceAllSimpleMeanDiff")
    run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    if (is_kk_design){
      inference_banner("SeqDesignInferenceAllKKCompoundMeanDiff")
      run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
      inference_banner("SeqDesignInferenceContinMultOLSKK")
      run_inference_checks(SeqDesignInferenceContinMultOLSKK$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
      inference_banner("SeqDesignInferenceContinMultGLS")
      run_inference_checks(SeqDesignInferenceContinMultGLS$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    } else {
      inference_banner("SeqDesignInferenceContinMultOLS")
      run_inference_checks(SeqDesignInferenceContinMultOLS$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    }
  }

  if (response_type == "incidence"){
    inference_banner("SeqDesignInferenceAllSimpleMeanDiff")
    run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    if (is_kk_design){
      inference_banner("SeqDesignInferenceAllKKCompoundMeanDiff")
      run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
      inference_banner("SeqDesignInferenceIncidKKClogit")
      run_inference_checks(SeqDesignInferenceIncidKKClogit$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    }
    inference_banner("SeqDesignInferenceIncidUnivLogRegr")
    run_inference_checks(SeqDesignInferenceIncidUnivLogRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    inference_banner("SeqDesignInferenceIncidMultiLogRegr")
    run_inference_checks(SeqDesignInferenceIncidMultiLogRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
  }

  if (response_type == "proportion"){
    inference_banner("SeqDesignInferenceAllSimpleMeanDiff")
    run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    if (is_kk_design){
      inference_banner("SeqDesignInferenceAllKKCompoundMeanDiff")
      run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    }    
    inference_banner("SeqDesignInferencePropUniBetaRegr")
    run_inference_checks(SeqDesignInferencePropUniBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    inference_banner("SeqDesignInferencePropMultiBetaRegr")
    run_inference_checks(SeqDesignInferencePropMultiBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
  }

  if (response_type == "count"){
    inference_banner("SeqDesignInferenceAllSimpleMeanDiff")
    run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    if (is_kk_design){
      inference_banner("SeqDesignInferenceAllKKCompoundMeanDiff")
      run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    }    
    inference_banner("SeqDesignInferenceCountUnivNegBinRegr")
    run_inference_checks(SeqDesignInferenceCountUnivNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    inference_banner("SeqDesignInferenceCountMultiNegBinRegr")
    run_inference_checks(SeqDesignInferenceCountMultiNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
  }

  if (response_type == "survival"){
    inference_banner("SeqDesignInferenceSurvivalRestrictedMeanDiff")
    run_inference_checks(SeqDesignInferenceSurvivalRestrictedMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    inference_banner("SeqDesignInferenceSurvivalKMDiff")
    run_inference_checks(SeqDesignInferenceSurvivalKMDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    inference_banner("SeqDesignInferenceSurvivalUniWeibullRegr")
    run_inference_checks(SeqDesignInferenceSurvivalUniWeibullRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    inference_banner("SeqDesignInferenceSurvivalMultiWeibullRegrv")
    run_inference_checks(SeqDesignInferenceSurvivalMultiWeibullRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    inference_banner("SeqDesignInferenceSurvivalUniCoxPHRegr")
    run_inference_checks(SeqDesignInferenceSurvivalUniCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
    inference_banner("SeqDesignInferenceSurvivalMultiCoxPHRegr")
    run_inference_checks(SeqDesignInferenceSurvivalMultiCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
  }
}


# response_type = "incidence"
# design_type = "CRD"
# dataset_name = "boston"

for (rep_curr in 1:Nrep) {
  message(paste0("\n\n====== RUNNING REPETITION ", rep_curr, " OF ", Nrep, " ======\n\n"))
  for (beta_T_iter_curr in seq_along(beta_T_values)){
    beta_T = beta_T_values[beta_T_iter_curr]
    message(paste0("\n\n=== Running tests with beta_T = [", beta_T, "] ===\n\n"))
    for (dataset_name in names(datasets_and_response_models)){
      message(paste0("\n\n=== Running tests for dataset: ", dataset_name, " ===\n\n"))
      for (response_type in c("continuous", "incidence", "proportion", "count", "survival")) {
        if (!(response_type %in% names(datasets_and_response_models[[dataset_name]]$y_original))) {
          message(paste0("\n\n  === Skipping response_type: ", response_type, " for dataset: ", dataset_name, " (not found in y_original) ===\n\n"))
          next
        }
        message(paste0("\n\n  === Running tests for response_type: ", response_type, " ===\n\n"))
        for (design_type in c("CRD", "iBCRD", "Efron", "KK14", "KK21", "KK21stepwise")) { #, "Atkinson"
          message(paste0("\n\n    === Running tests for design: ", design_type, " ==="))
          run_tests_for_response(response_type, design_type = design_type, dataset_name = dataset_name)
          message(paste0("\n\n  === Finished tests for design_type: ", design_type, " ===\n\n"))
        }
        message(paste0("\n\n  === Finished tests for response_type: ", response_type, " ===\n\n"))
      }
      message(paste0("\n\n  === Finished tests for dataset: ", dataset_name, " ===\n\n"))
    }
  }
}
message("\n\n----------------------All tests complete!")
write_results_if_needed(force = TRUE)

rm(list=ls()) #; rstudioapi::restartSession()
