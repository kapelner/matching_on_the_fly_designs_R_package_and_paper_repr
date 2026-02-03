#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.5/ SeqExpMatch/
rm(list = ls())
set.seed(1)
pacman::p_load(SeqExpMatch, doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
max_n_dataset = 200
source("package_tests/_dataset_load.R")
# options(warn=2)

NUM_CORES = 1 # Force serial execution for debugging
prob_censoring = 0.3
nsim_exact_test = 101
pval_epsilon = 0.01
test_compute_confidence_interval_rand = TRUE #too slow for now

results_file = "package_tests/simple_tests_results.csv"
results_dt = data.table(
  dataset = character(),
  response_type = character(),
  design_type = character(),
  function_run = character(),
  result = character()
)
results_write_every = 50L
results_last_written = 0L

write_results_if_needed = function(force = FALSE){
  if (force || (nrow(results_dt) - results_last_written) >= results_write_every){
    data.table::fwrite(results_dt, results_file)
    results_last_written <<- nrow(results_dt)
  }
}

record_result = function(dataset_name, response_type, design_type, function_run, result){
  result_str = if (is.null(result)) {
    NA_character_
  } else if (length(result) == 0) {
    ""
  } else if (length(result) == 1) {
    as.character(result)
  } else {
    paste(as.character(result), collapse = " ")
  }
  results_dt <<- data.table::rbindlist(list(
    results_dt,
    data.table(
      dataset = dataset_name,
      response_type = response_type,
      design_type = design_type,
      function_run = function_run,
      result = result_str
    )
  ))
  write_results_if_needed()
}

run_inference_checks = function(seq_des_inf, response_type, design_type, dataset_name){
  safe_call = function(label, expr){
    tryCatch({
        result <- expr
        cat("  ", result, "\n")
        record_result(dataset_name, response_type, design_type, label, result)
        result
      }, error = function(e){
      message("Skipping ", label, ": ", e$message)
      record_result(dataset_name, response_type, design_type, label, NA_character_)
      NULL
    })
  }

  message("    Calling compute_treatment_estimate()")
  safe_call("compute_treatment_estimate",
    seq_des_inf$compute_treatment_estimate())
  message("    Calling compute_mle_two_sided_pval_for_treatment_effect()")
  safe_call("compute_mle_two_sided_pval_for_treatment_effect",
            seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect())
  message("    Calling compute_mle_confidence_interval()")
  safe_call("compute_mle_confidence_interval",
            seq_des_inf$compute_mle_confidence_interval(0.05))
  message("    Calling compute_bootstrap_confidence_interval()")
  safe_call("compute_bootstrap_confidence_interval",
            seq_des_inf$compute_bootstrap_confidence_interval(B = nsim_exact_test, na.rm = TRUE))
  message("    Calling compute_bootstrap_two_sided_pval()")
  safe_call("compute_bootstrap_two_sided_pval",
            seq_des_inf$compute_bootstrap_two_sided_pval(B = nsim_exact_test, na.rm = TRUE))
  if (response_type %in% c("continuous", "survival")){
    message("    Calling compute_two_sided_pval_for_treatment_effect_rand()")
    safe_call("compute_two_sided_pval_for_treatment_effect_rand",
            seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = nsim_exact_test, show_progress = FALSE))
    message("    Calling compute_two_sided_pval_for_treatment_effect_rand(delta=0.5)")
    is_cox = inherits(seq_des_inf, c("SeqDesignInferenceSurvivalUniCoxPHRegr", "SeqDesignInferenceSurvivalMultiCoxPHRegr"))
     safe_call("compute_two_sided_pval_for_treatment_effect_rand(delta=0.5)",
            seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(
              nsim_exact_test = nsim_exact_test,
              delta = 0.5,
              transform_responses = ifelse(response_type == "survival" && !is_cox, "log", "none"),
              show_progress = FALSE
            ))
  }

  if (test_compute_confidence_interval_rand & response_type %in% c("continuous")){ #,  "proportion", "survival"
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
  message("    Calling compute_two_sided_pval_for_treatment_effect_rand(custom)")
  safe_call("compute_two_sided_pval_for_treatment_effect_rand(custom)",
            seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = nsim_exact_test, show_progress = FALSE))
  if (test_compute_confidence_interval_rand & response_type %in% c("continuous")){ #, "proportion", "survival"
    message("    Calling compute_confidence_interval_rand(custom)")
    safe_call("compute_confidence_interval_rand(custom)",
              seq_des_inf$compute_confidence_interval_rand(nsim_exact_test = nsim_exact_test, pval_epsilon = pval_epsilon, show_progress = FALSE))
  }
  seq_des_inf$set_custom_randomization_statistic_function(NULL)
}

run_tests_for_response = function(response_type, design_type, dataset_name){
  D = datasets_and_response_models[[dataset_name]]
  n = nrow(D$X)
  dead = rep(1, n)
  y = D$y_original[[response_type]]

  if (response_type == "survival"){
    for (i in 1 : n){
      if (runif(1) < prob_censoring){
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
    seq_des_obj$add_subject_to_experiment_and_assign(D$X[t, ])
    seq_des_obj$add_subject_response(t, y[t], dead[t])
  }

  is_kk_design = design_type %in% c("KK21", "KK21stepwise", "KK14")
  if (response_type == "continuous"){
    message("\n\n    == Inference: SeqDesignInferenceAllSimpleMeanDiff\n")
    run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    if (is_kk_design){
      message("\n\n    == Inference: SeqDesignInferenceAllKKCompoundMeanDiff\n")
      run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
      message("\n\n    == Inference: SeqDesignInferenceContinMultOLSKK\n")
      run_inference_checks(SeqDesignInferenceContinMultOLSKK$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    } else {
      message("\n\n    == Inference: SeqDesignInferenceContinMultOLS\n")
      run_inference_checks(SeqDesignInferenceContinMultOLS$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    }
  }

  if (response_type == "incidence"){
    message("\n\n  == Inference: SeqDesignInferenceAllSimpleMeanDiff\n")
    run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    if (is_kk_design){
      message("\n\n    == Inference: SeqDesignInferenceAllKKCompoundMeanDiff\n")
      run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    }    
    message("\n\n  == Inference: SeqDesignInferenceIncidUnivLogRegr\n")
    run_inference_checks(SeqDesignInferenceIncidUnivLogRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    message("\n\n  == Inference: SeqDesignInferenceIncidMultiLogRegr\n")
    run_inference_checks(SeqDesignInferenceIncidMultiLogRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
  }

  if (response_type == "proportion"){
    message("\n\n  == Inference: SeqDesignInferenceAllSimpleMeanDiff\n")
    run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    if (is_kk_design){
      message("\n\n    == Inference: SeqDesignInferenceAllKKCompoundMeanDiff\n")
      run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    }    
    message("\n\n  == Inference: SeqDesignInferencePropUniBetaRegr\n")
    run_inference_checks(SeqDesignInferencePropUniBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    message("\n\n  == Inference: SeqDesignInferencePropMultiBetaRegr\n")
    run_inference_checks(SeqDesignInferencePropMultiBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
  }

  if (response_type == "count"){
    message("\n\n  == Inference: SeqDesignInferenceAllSimpleMeanDiff\n")
    run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    if (is_kk_design){
      message("\n\n    == Inference: SeqDesignInferenceAllKKCompoundMeanDiff\n")
      run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    }    
    message("\n\n  == Inference: SeqDesignInferenceCountUnivNegBinRegr\n")
    run_inference_checks(SeqDesignInferenceCountUnivNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    message("\n\n  == Inference: SeqDesignInferenceCountMultiNegBinRegr\n")
    run_inference_checks(SeqDesignInferenceCountMultiNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
  }

  if (response_type == "survival"){
    message("\n\n  == Inference: SeqDesignInferenceSurvivalRestrictedMeanDiff\n")
    run_inference_checks(SeqDesignInferenceSurvivalRestrictedMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    message("\n\n  == Inference: SeqDesignInferenceSurvivalKMDiff\n")
    run_inference_checks(SeqDesignInferenceSurvivalKMDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    message("\n\n  == Inference: SeqDesignInferenceSurvivalUniWeibullRegr\n")
    run_inference_checks(SeqDesignInferenceSurvivalUniWeibullRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    message("\n\n  == Inference: SeqDesignInferenceSurvivalMultiWeibullRegrv")
    run_inference_checks(SeqDesignInferenceSurvivalMultiWeibullRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    message("\n\n  == Inference: SeqDesignInferenceSurvivalUniCoxPHRegr\n")
    run_inference_checks(SeqDesignInferenceSurvivalUniCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
    message("\n\n  == Inference: SeqDesignInferenceSurvivalMultiCoxPHRegr\n")
    run_inference_checks(SeqDesignInferenceSurvivalMultiCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name)
  }
}


# response_type = "incidence"
# design_type = "CRD"
# dataset_name = "boston"

for (dataset_name in names(datasets_and_response_models)){
  message(paste0("\n\n=== Running tests for dataset: ", dataset_name, " ===\n\n"))
  for (response_type in c("continuous", "incidence", "proportion", "count", "survival")) {
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
message("\n\n----------------------All tests complete!")
write_results_if_needed(force = TRUE)

rm(list=ls()); rstudioapi::restartSession()
