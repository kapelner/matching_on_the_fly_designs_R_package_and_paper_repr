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
  betaT = 1
  sd_noise = 0.1
  
  for (i in 1:n) {
    w_i = des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
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
  
  for (inf_info in inference_classes) {
    inf_cls = if (is.list(inf_info)) inf_info[[1]] else inf_info
    inf_args = if (is.list(inf_info) && length(inf_info) > 1) inf_info[-1] else list()
    
    cls_name = if (is.character(inf_cls)) inf_cls else inf_cls$classname
    cat("\n--- Inference:", cls_name, if (length(inf_args)) paste0("(", paste(names(inf_args), inf_args, sep="=", collapse=", "), ")") else "", "---\n")
    
    inf = do.call(inf_cls$new, c(list(des_obj = des), inf_args))
    
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
  list(InferenceContinOLS, include_covariates = FALSE),
  list(InferenceContinOLS, include_covariates = TRUE),
  list(InferenceContinLin, include_covariates = TRUE),
  list(InferenceContinQuantileRegr, include_covariates = TRUE),
  list(InferenceContinRobustRegr, include_covariates = TRUE),
  list(InferenceContinKKRobustRegrIVWC, include_covariates = TRUE),
  list(InferenceContinKKRobustRegrOneLik, include_covariates = TRUE),
  list(InferenceContinKKGLMM, include_covariates = TRUE)
))

##### response_type = incidence
run_tests_for_response("incidence", list(
  list(InferenceIncidLogRegr, include_covariates = FALSE),
  list(InferenceIncidLogRegr, include_covariates = TRUE),
  list(InferenceIncidLogBinomial, include_covariates = TRUE),
  list(InferenceIncidModifiedPoisson, include_covariates = TRUE),
  list(InferenceIncidRiskDiff, include_covariates = TRUE),
  list(InferenceIncidBinomialIdentityRiskDiff, include_covariates = TRUE),
  list(InferenceIncidGCompRiskDiff, include_covariates = TRUE),
  list(InferenceIncidGCompRiskRatio, include_covariates = TRUE),
  list(InferenceIncidKKClogitIVWC, include_covariates = TRUE),
  list(InferenceIncidKKClogitOneLik, include_covariates = TRUE),
  list(InferenceIncidKKGEE, include_covariates = TRUE),
  list(InferenceIncidKKGLMM, include_covariates = TRUE)
))

##### response_type = proportion
run_tests_for_response("proportion", list(
  InferenceAllSimpleMeanDiff,
  list(InferencePropBetaRegr, include_covariates = TRUE),
  list(InferencePropFractionalLogit, include_covariates = TRUE),
  list(InferencePropGCompMeanDiff, include_covariates = TRUE),
  list(InferencePropKKGEE, include_covariates = TRUE),
  list(InferencePropKKGLMM, include_covariates = TRUE)
))

##### response_type = count
run_tests_for_response("count", list(
  list(InferenceCountPoisson, include_covariates = TRUE),
  list(InferenceCountNegBin, include_covariates = TRUE),
  list(InferenceCountKKCPoissonIVWC, include_covariates = TRUE)
))

##### response_type = survival
run_tests_for_response("survival", list(
  InferenceSurvivalLogRank,
  list(InferenceSurvivalCoxPHRegr, include_covariates = TRUE),
  list(InferenceSurvivalKKLWACoxIVWC, include_covariates = TRUE)
))

##### response_type = ordinal
run_tests_for_response("ordinal", list(
  InferenceOrdinalJonckheereTerpstraTest
))
