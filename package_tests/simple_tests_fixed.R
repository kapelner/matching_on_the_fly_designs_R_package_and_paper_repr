#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.5/ SeqExpMatch/
rm(list = ls())
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
.libPaths(c("Rlib", .libPaths()))
library(SeqExpMatch)
debug_mode = nzchar(Sys.getenv("SEQEXPMATCH_DEBUG"))
if (debug_mode){
  options(error = function(){
    traceback(2)
    q(status = 1)
  })
} else {
  options(error = NULL)
}
# options(warn=2)
set.seed(1)

max_n_dataset = 500
source("_dataset_load.R")


D = datasets_and_response_models$boston

#try to create a CRD design
n = min(nrow(D$X), 200)
dead = rep(1, n)
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

run_tests_for_response = function(response_type){
  if (debug_mode){
    message("Running response type: ", response_type)
  }
  y = D$y_original[[response_type]]
  dead = rep(1, n)

  if (response_type == "survival"){
    prob_censoring = 0.3
    for (i in 1 : n){
      if (runif(1) < prob_censoring){
        y[i] = runif(1, 0, y[i])
        dead[i] = 0
      }
    }
  }

  seq_des_obj = SeqDesignKK21$new(response_type = response_type, n = n)
  for (t in 1 : n){
    seq_des_obj$add_subject_to_experiment_and_assign(D$X[t, ])
    seq_des_obj$add_subject_response(t, y[t], dead[t])
  }

  if (response_type == "continuous"){
    if (debug_mode){
      message("  Inference: SeqDesignInferenceAllSimpleMeanDiff")
    }
    run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj), response_type)
    if (debug_mode){
      message("  Inference: SeqDesignInferenceAllKKCompoundMeanDiff")
    }
    run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj), response_type)
    if (debug_mode){
      message("  Inference: SeqDesignInferenceContinMultOLSKK")
    }
    run_inference_checks(SeqDesignInferenceContinMultOLSKK$new(seq_des_obj), response_type)
  }

  if (response_type == "incidence"){
    if (debug_mode){
      message("  Inference: SeqDesignInferenceIncidUnivLogRegr")
    }
    run_inference_checks(SeqDesignInferenceIncidUnivLogRegr$new(seq_des_obj), response_type)
    if (debug_mode){
      message("  Inference: SeqDesignInferenceIncidMultiLogRegr")
    }
    run_inference_checks(SeqDesignInferenceIncidMultiLogRegr$new(seq_des_obj), response_type)
  }

  if (response_type == "proportion"){
    if (debug_mode){
      message("  Inference: SeqDesignInferencePropUniBetaRegr")
    }
    run_inference_checks(SeqDesignInferencePropUniBetaRegr$new(seq_des_obj), response_type)
    if (debug_mode){
      message("  Inference: SeqDesignInferencePropMultiBetaRegr")
    }
    run_inference_checks(SeqDesignInferencePropMultiBetaRegr$new(seq_des_obj), response_type)
  }

  if (response_type == "count"){
    if (debug_mode){
      message("  Inference: SeqDesignInferenceCountUnivNegBinRegr")
    }
    run_inference_checks(SeqDesignInferenceCountUnivNegBinRegr$new(seq_des_obj), response_type)
    if (debug_mode){
      message("  Inference: SeqDesignInferenceCountMultiNegBinRegr")
    }
    run_inference_checks(SeqDesignInferenceCountMultiNegBinRegr$new(seq_des_obj), response_type)
  }

  if (response_type == "survival"){
    if (debug_mode){
      message("  Inference: SeqDesignInferenceSurvivalRestrictedMeanDiff")
    }
    run_inference_checks(SeqDesignInferenceSurvivalRestrictedMeanDiff$new(seq_des_obj), response_type)
    if (debug_mode){
      message("  Inference: SeqDesignInferenceSurvivalKMDiff")
    }
    run_inference_checks(SeqDesignInferenceSurvivalKMDiff$new(seq_des_obj), response_type)
    if (debug_mode){
      message("  Inference: SeqDesignInferenceSurvivalUniWeibullRegr")
    }
    run_inference_checks(SeqDesignInferenceSurvivalUniWeibullRegr$new(seq_des_obj), response_type)
    if (debug_mode){
      message("  Inference: SeqDesignInferenceSurvivalMultiWeibullRegr")
    }
    run_inference_checks(SeqDesignInferenceSurvivalMultiWeibullRegr$new(seq_des_obj), response_type)
    if (debug_mode){
      message("  Inference: SeqDesignInferenceSurvivalUniCoxPHRegr")
    }
    run_inference_checks(SeqDesignInferenceSurvivalUniCoxPHRegr$new(seq_des_obj), response_type)
    if (debug_mode){
      message("  Inference: SeqDesignInferenceSurvivalMultiCoxPHRegr")
    }
    run_inference_checks(SeqDesignInferenceSurvivalMultiCoxPHRegr$new(seq_des_obj), response_type)
  }
}

run_tests_for_response("continuous")
run_tests_for_response("incidence")
run_tests_for_response("proportion")
run_tests_for_response("count")
run_tests_for_response("survival")
