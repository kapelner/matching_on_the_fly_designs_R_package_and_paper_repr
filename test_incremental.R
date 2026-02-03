#!/usr/bin/env Rscript
# This will test the exact lines from simple_tests.R one by one

cat("Step 1: rm(list = ls())\n")
rm(list = ls())

cat("Step 2: Loading packages\n")
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)

cat("Step 3: Script path detection\n")
file_arg = "--file="
script_path = sub(file_arg, "", commandArgs(trailingOnly = FALSE)[grep(file_arg, commandArgs(trailingOnly = FALSE))])
test_dir = dirname(normalizePath(script_path))
lib_path = normalizePath(file.path(test_dir, "..", "Rlib"), mustWork = FALSE)
if (dir.exists(lib_path)){
  .libPaths(c(lib_path, .libPaths()))
}

cat("Step 4: Loading SeqExpMatch\n")
library(SeqExpMatch)

cat("Step 5: Setting debug_mode\n")
debug_mode = nzchar(Sys.getenv("SEQEXPMATCH_DEBUG"))
if (debug_mode){
  options(error = function(){
    traceback(2)
    q(status = 1)
  })
} else {
  options(error = NULL)
}

cat("Step 6: set.seed(1)\n")
set.seed(1)

cat("Step 7: max_n_dataset = 500\n")
max_n_dataset = 500

cat("Step 8: source('package_tests/_dataset_load.R')\n")
source("package_tests/_dataset_load.R")

cat("Step 9: D = datasets_and_response_models$boston\n")
D = datasets_and_response_models$boston

cat("Step 10: Setting n and dead\n")
n = min(nrow(D$X), 200)
dead = rep(1, n)

cat("Step 11: Defining run_inference_checks function\n")
run_inference_checks = function(seq_des_inf, response_type){
  safe_call = function(label, expr){
    tryCatch(expr, error = function(e){
      message("Skipping ", label, ": ", e$message)
      NULL
    })
  }

  seq_des_inf$compute_treatment_estimate()
  safe_call("compute_mle_two_sided_pval_for_treatment_effect",
            seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect())
  safe_call("compute_mle_confidence_interval",
            seq_des_inf$compute_mle_confidence_interval(0.05))
  safe_call("compute_bootstrap_confidence_interval",
            seq_des_inf$compute_bootstrap_confidence_interval(B = 50, na.rm = TRUE))
  safe_call("compute_bootstrap_two_sided_pval",
            seq_des_inf$compute_bootstrap_two_sided_pval(B = 50, na.rm = TRUE))
  safe_call("compute_two_sided_pval_for_treatment_effect_rand",
            seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = 50))
  if (response_type %in% c("continuous", "proportion", "survival")){
     safe_call("compute_two_sided_pval_for_treatment_effect_rand(delta=0.5)",
            seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = 50, delta = 0.5))
  }
  safe_call("compute_two_sided_pval_for_treatment_effect_rand(na.rm=TRUE)",
            seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = 50, na.rm = TRUE))

  if (response_type %in% c("continuous", "proportion", "survival")){
    safe_call("compute_confidence_interval_rand",
              seq_des_inf$compute_confidence_interval_rand(nsim_exact_test = 50, pval_epsilon = 0.05))
    safe_call("compute_confidence_interval_rand(alpha=0.01)",
              seq_des_inf$compute_confidence_interval_rand(nsim_exact_test = 50, pval_epsilon = 0.05, alpha = 0.01))
  }
  seq_des_inf$set_custom_randomization_statistic_function(function(){
    yTs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 1]
    yCs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 0]
    (mean(yTs) - mean(yCs)) / sqrt(var(yTs) / length(yTs) + var(yCs) / length(yCs))
  })
  safe_call("compute_two_sided_pval_for_treatment_effect_rand(custom)",
            seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = 50))
  if (response_type %in% c("continuous", "proportion", "survival")){
    safe_call("compute_confidence_interval_rand(custom)",
              seq_des_inf$compute_confidence_interval_rand(nsim_exact_test = 50, pval_epsilon = 0.05))
  }
  seq_des_inf$set_custom_randomization_statistic_function(NULL)
}

cat("Step 12: Defining run_tests_for_response function\n")
run_tests_for_response = function(response_type){
  cat("  Inside run_tests_for_response, response_type:", response_type, "\n")
  y = D$y_original[[response_type]]
  cat("  Got y, length:", length(y), "\n")
  dead = rep(1, n)
  cat("  Got dead\n")

  if (response_type == "survival"){
    prob_censoring = 0.3
    for (i in 1 : n){
      if (runif(1) < prob_censoring){
        y[i] = runif(1, 0, y[i])
        dead[i] = 0
      }
    }
  }

  cat("  Creating SeqDesignKK21 object\n")
  seq_des_obj = SeqDesignKK21$new(response_type = response_type, n = n)
  cat("  Created SeqDesignKK21 object\n")

  cat("  Adding subjects...\n")
  for (t in 1 : n){
    if (t %% 50 == 0) cat("    Added", t, "subjects\n")
    seq_des_obj$add_subject_to_experiment_and_assign(D$X[t, ])
    seq_des_obj$add_subject_response(t, y[t], dead[t])
  }
  cat("  Added all", n, "subjects\n")

  if (response_type == "continuous"){
    cat("  About to run inference checks\n")
    cat("  Creating SeqDesignInferenceAllSimpleMeanDiff\n")
    run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj), response_type)
    cat("  Done with SeqDesignInferenceAllSimpleMeanDiff\n")

    cat("  Creating SeqDesignInferenceAllKKCompoundMeanDiff\n")
    run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj), response_type)
    cat("  Done with SeqDesignInferenceAllKKCompoundMeanDiff\n")

    cat("  Creating SeqDesignInferenceContinMultOLSKK\n")
    run_inference_checks(SeqDesignInferenceContinMultOLSKK$new(seq_des_obj), response_type)
    cat("  Done with SeqDesignInferenceContinMultOLSKK\n")
  }

  cat("  Finished run_tests_for_response\n")
}

cat("Step 13: Running run_tests_for_response('continuous')\n")
run_tests_for_response("continuous")

cat("SUCCESS: All steps completed!\n")
