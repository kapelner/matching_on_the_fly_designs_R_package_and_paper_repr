library(EDI)

# Settings
n = 100
p = 5
data_type = "linear"
design_cls = DesignSeqOneByOneKK21
r = 50 # number of bootstrap/randomization iterations for speed

# 1. Generate Data
dataset = generate_covariate_dataset(n, p, data_type)
X = dataset$X
y_cont = dataset$y_cont

run_tests_for_response = function(response_type, inference_classes) {
  cat("\n\n###########################################################################\n")
  cat("##### response_type =", response_type, "\n")
  cat("###########################################################################\n")
  
  y_base = transform_cont_y_based_on_response_type(y_cont, response_type)
  
  # Initialize Design
  des = design_cls$new(response_type = response_type, n = n)
  
  # Fill Design (Sequential)
  # For simulation purposes, we need to apply treatment effect and noise like SimulationFramework does
  betaT = 1
  sd_noise = 0.1
  
  for (i in 1:n) {
    w_i = des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
    
    # Simple version of .apply_one logic
    bt = if (w_i == 1) betaT else 0
    eps = rnorm(1, 0, sd_noise)
    
    y_i = switch(response_type,
      continuous = y_base[i] + bt + eps,
      incidence  = {
        p_b = if (is.finite(y_base[i]) && y_base[i] >= 0 && y_base[i] <= 1) y_base[i] else plogis(y_base[i])
        p_b = pmin(0.95, pmax(0.05, p_b))
        as.numeric(rbinom(1, 1, plogis(qlogis(p_b) + bt + eps)))
      },
      proportion = pmin(1 - 1e-9, pmax(1e-9, y_base[i] + bt + eps)),
      count = {
        lam = pmax(.Machine$double.eps, y_base[i] * exp(bt + eps))
        as.numeric(rpois(1, lambda = lam))
      },
      survival = pmax(.Machine$double.eps, y_base[i] * exp(bt + eps)),
      ordinal  = as.integer(max(1, round(y_base[i] + bt + eps))),
      stop("Unknown response_type")
    )
    
    des$add_one_subject_response(i, y_i, dead = 1)
  }
  
  for (inf_cls in inference_classes) {
    cat("\n--- Inference:", inf_cls$classname, "---\n")
    inf = inf_cls$new(des)
    
    # 1. Treatment Estimate
    tryCatch({
      est = inf$compute_treatment_estimate()
      cat("  Estimate:", est, "\n")
    }, error = function(e) cat("  Estimate Error:", e$message, "\n"))
    
    # 2. Asymptotic Inference
    if (inherits(inf, "InferenceAsymp")) {
      tryCatch({
        cat("  Asymp P-val:", inf$compute_asymp_two_sided_pval_for_treatment_effect(), "\n")
        cat("  Asymp CI:", paste(inf$compute_asymp_confidence_interval(), collapse = ", "), "\n")
      }, error = function(e) cat("  Asymp Error:", e$message, "\n"))
    }
    
    # 3. Bootstrap Inference
    if (inherits(inf, "InferenceBoot")) {
      tryCatch({
        cat("  Boot P-val:", inf$compute_bootstrap_two_sided_pval(B = r), "\n")
        cat("  Boot CI:", paste(inf$compute_bootstrap_confidence_interval(B = r), collapse = ", "), "\n")
      }, error = function(e) cat("  Boot Error:", e$message, "\n"))
    }
    
    # 4. Randomization Inference
    if (inherits(inf, "InferenceRand")) {
      tryCatch({
        cat("  Rand P-val:", inf$compute_two_sided_pval_for_treatment_effect_rand(r = r), "\n")
        if (inherits(inf, "InferenceRandCI")) {
          cat("  Rand CI:", paste(inf$compute_confidence_interval_rand(r = r), collapse = ", "), "\n")
        }
      }, error = function(e) cat("  Rand Error:", e$message, "\n"))
    }
  }
}

##### response_type = continuous
run_tests_for_response("continuous", list(
  InferenceAllSimpleMeanDiff,
  InferenceContinMultOLS,
  InferenceContinMultOLSKKIVWC
))

##### response_type = incidence
run_tests_for_response("incidence", list(
  InferenceIncidUnivLogRegr,
  InferenceIncidMultiLogRegr,
  InferenceIncidUnivKKClogitIVWC
))

##### response_type = proportion
run_tests_for_response("proportion", list(
  InferenceAllSimpleMeanDiff,
  InferencePropUniBetaRegr,
  InferencePropMultiKKGEE
))

##### response_type = count
run_tests_for_response("count", list(
  InferenceCountUnivPoissonRegr,
  InferenceCountMultiNegBinRegr,
  InferenceCountPoissonUnivKKCPoissonIVWC
))

##### response_type = survival
run_tests_for_response("survival", list(
  InferenceSurvivalLogRank,
  InferenceSurvivalUniCoxPHRegr,
  InferenceSurvivalUnivKKLWACoxIVWC
))

##### response_type = ordinal
run_tests_for_response("ordinal", list(
  InferenceOrdinalUniPropOddsRegr,
  InferenceOrdinalJonckheereTerpstraTest,
  InferenceOrdinalUnivKKGEE
))
