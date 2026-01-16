if (!requireNamespace("pkgload", quietly = TRUE)) {
  install.packages("pkgload", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(survival))

pkgload::load_all("SeqExpMatch")
args <- commandArgs(trailingOnly = TRUE)
verbose <- !("--quiet" %in% args)
vbcat <- function(...) {
  #if (verbose) {
    cat(...)
  #}
}
vbcat("Package loaded from:", path.package("SeqExpMatch"), "\n")

# Helpers to safely run inference methods and normalize confidence interval outputs
safe_run_metric <- function(expr) {
  tryCatch(
    list(value = expr(), error = NA_character_),
    error = function(e) list(value = NA, error = conditionMessage(e))
  )
}

extract_ci_bounds <- function(value) {
  if (is.numeric(value) && length(value) >= 2) {
    return(list(lower = as.numeric(value[1]), upper = as.numeric(value[2])))
  }
  list(lower = NA_real_, upper = NA_real_)
}

# Set seed for reproducibility
set.seed(123)

designs <- c("SeqDesignCRD", "SeqDesigniBCRD", "SeqDesignEfron", "SeqDesignAtkinson", "SeqDesignKK14", "SeqDesignKK21", "SeqDesignKK21stepwise")
response_types <- c("continuous", "incidence", "count", "proportion", "survival")

inf_by_res_type <- list(
  continuous = c(
    "SeqDesignInferenceAllSimpleMeanDiff",
    "SeqDesignInferenceAllKKCompoundMeanDiff",
    "SeqDesignInferenceContinMultOLS",
    "SeqDesignInferenceContinMultOLSKK",
    "SeqDesignInferenceBaiAdjustedTKK14",
    "SeqDesignInferenceBaiAdjustedTKK21"
  ),
  incidence = c(
    "SeqDesignInferenceAllSimpleMeanDiff",
    "SeqDesignInferenceAllKKCompoundMeanDiff",
    "SeqDesignInferenceIncidUnivLogRegr",
    "SeqDesignInferenceIncidMultiLogRegr"
  ),
  count = c(
    "SeqDesignInferenceAllSimpleMeanDiff",
    "SeqDesignInferenceAllKKCompoundMeanDiff",
    "SeqDesignInferenceCountUnivNegBinRegr",
    "SeqDesignInferenceCountMultiNegBinRegr"
  ),
  proportion = c(
    "SeqDesignInferenceAllSimpleMeanDiff",
    "SeqDesignInferenceAllKKCompoundMeanDiff",
    "SeqDesignInferencePropUniBetaRegr",
    "SeqDesignInferencePropMultiBetaRegr"
  ),
  survival = c(
    "SeqDesignInferenceSurvivalKMDiff",
    "SeqDesignInferenceSurvivalRestrictedMeanDiff",
    "SeqDesignInferenceSurvivalUniWeibullRegr",
    "SeqDesignInferenceSurvivalMultiWeibullRegr",
    "SeqDesignInferenceSurvivalUniCoxPHRegr",
    "SeqDesignInferenceSurvivalMultiCoxPHRegr"
  )
)

inference_design_requirements <- list(
  SeqDesignInferenceBaiAdjustedTKK14 = "SeqDesignKK14",
  SeqDesignInferenceBaiAdjustedTKK21 = c("SeqDesignKK21", "SeqDesignKK21stepwise")
)

n <- 100
p <- 5
beta_T = 1

# Generate synthetic covariate data
X_mat <- matrix(rnorm(n * p), nrow = n)
colnames(X_mat) <- paste0("x", 1:p)
X_dt <- as.data.table(X_mat)

results <- list()

for (rt in response_types) {
  vbcat("\nProcessing response type:", rt, "\n")
  
  # Generate responses
  y_controls <- rnorm(n)
  deads <- if (rt == "survival") rbinom(n, 1, 0.8) else rep(1, n)
  
  for (des_name in designs) {
    vbcat("  Testing design:", des_name, "\n")
    
    # Initialize design
    # Some designs might have different constructor arguments, but most follow SeqDesign signature
    des_class <- get(des_name)
    des_obj <- des_class$new(n = n, response_type = rt, verbose = verbose)
    
    # Run design
    for (t in 1:n) {
      w_t = des_obj$add_subject_to_experiment_and_assign(X_dt[t, ])
      eta_t = y_controls[t] + beta_T * w_t
      y_t <- switch(rt,
                continuous = rnorm(1, eta_t, 1),
                incidence = rbinom(1, 1, inv_logit(eta_t)),
                count = rpois(1, exp(eta_t)),
                proportion = inv_logit(rnorm(1, eta_t, 1)),
                survival = rexp(1, exp(eta_t))
              )     
      des_obj$add_subject_response(t, y_t, dead = deads[t])
    }
    
    # Test each inference method
    inf_methods <- inf_by_res_type[[rt]]
    for (inf_name in inf_methods) {
      vbcat("    Inference method:", inf_name, "... ")
      req_class <- inference_design_requirements[[inf_name]]
      skip_reason <- NULL
      if (!is.null(req_class)) {
        allowed_classes <- if (is.character(req_class)) req_class else as.character(req_class)
        if (!any(vapply(allowed_classes, function(cl) inherits(des_obj, cl), logical(1)))) {
          skip_reason <- paste("Skipped: requires one of", paste(allowed_classes, collapse = ", "))
        }
      }
      if (is.null(skip_reason) &&
          inf_name == "SeqDesignInferenceContinMultOLSKK" &&
          !any(vapply(c("SeqDesignKK14", "SeqDesignKK21"), function(cl) inherits(des_obj, cl), logical(1)))) {
        skip_reason <- "Skipped: design must inherit SeqDesignKK14 or SeqDesignKK21"
      }
      if (!is.null(skip_reason)) {
        vbcat(skip_reason, "\n")
        results[[length(results) + 1]] <- data.table(
          response_type = rt,
          design = des_name,
          inference = inf_name,
          status = skip_reason,
          time_ms = NA_real_,
          estimate = NA_real_,
          pval = NA_real_,
          bootstrap_ci_lower = NA_real_,
          bootstrap_ci_upper = NA_real_,
          bootstrap_ci_error = NA_character_,
          bootstrap_pval = NA_real_,
          bootstrap_pval_error = NA_character_,
          rand_pval = NA_real_,
          rand_pval_error = NA_character_,
          rand_ci_lower = NA_real_,
          rand_ci_upper = NA_real_,
          rand_ci_error = NA_character_,
          mle_ci_lower = NA_real_,
          mle_ci_upper = NA_real_,
          mle_ci_error = NA_character_,
          mle_pval = NA_real_,
          mle_pval_error = NA_character_
        )
        next
      }
      inf_class <- get(inf_name)

      # Use tryCatch because some inference methods have strict design requirements
      inf_start_time <- Sys.time()
      inf_obj_res <- tryCatch({
        inf_obj <- inf_class$new(des_obj, verbose = verbose)

        # Benchmark estimate computation
        est <- inf_obj$compute_treatment_estimate()

        bootstrap_ci_res <- safe_run_metric(function() inf_obj$compute_bootstrap_confidence_interval())
        bootstrap_ci_bounds <- extract_ci_bounds(bootstrap_ci_res$value)
        bootstrap_pval_res <- safe_run_metric(function() inf_obj$compute_bootstrap_two_sided_pval())

        rand_pval_res <- safe_run_metric(function() inf_obj$compute_two_sided_pval_for_treatment_effect_rand())
        rand_ci_res <- safe_run_metric(function() inf_obj$compute_confidence_interval_rand())
        rand_ci_bounds <- extract_ci_bounds(rand_ci_res$value)

        if ("compute_mle_confidence_interval" %in% names(inf_obj)) {
          mle_ci_res <- safe_run_metric(function() inf_obj$compute_mle_confidence_interval())
        } else {
          mle_ci_res <- list(value = NA, error = NA_character_)
        }
        mle_ci_bounds <- extract_ci_bounds(mle_ci_res$value)

        if ("compute_mle_two_sided_pval_for_treatment_effect" %in% names(inf_obj)) {
          mle_pval_res <- safe_run_metric(function() inf_obj$compute_mle_two_sided_pval_for_treatment_effect())
        } else {
          mle_pval_res <- list(value = NA, error = NA_character_)
        }

        list(
          est = est,
          bootstrap_ci_lower = bootstrap_ci_bounds$lower,
          bootstrap_ci_upper = bootstrap_ci_bounds$upper,
          bootstrap_ci_error = bootstrap_ci_res$error,
          bootstrap_pval = bootstrap_pval_res$value,
          bootstrap_pval_error = bootstrap_pval_res$error,
          rand_pval = rand_pval_res$value,
          rand_pval_error = rand_pval_res$error,
          rand_ci_lower = rand_ci_bounds$lower,
          rand_ci_upper = rand_ci_bounds$upper,
          rand_ci_error = rand_ci_res$error,
          mle_ci_lower = mle_ci_bounds$lower,
          mle_ci_upper = mle_ci_bounds$upper,
          mle_ci_error = mle_ci_res$error,
          mle_pval = mle_pval_res$value,
          mle_pval_error = mle_pval_res$error,
          pval = mle_pval_res$value,
          status = "Success"
        )
      }, error = function(e) {
        # cat("\nTraceback for error:", e$message, "\n")
        # traceback(e)
        list(status = paste("Error:", e$message))
      })
      inf_end_time <- Sys.time()
      
      vbcat(inf_obj_res$status, "\n")
      
      results[[length(results) + 1]] <- data.table(
        response_type = rt,
        design = des_name,
        inference = inf_name,
        status = inf_obj_res$status,
        time_ms = as.numeric(difftime(inf_end_time, inf_start_time, units = "secs")) * 1000,
        estimate = if (inf_obj_res$status == "Success") inf_obj_res$est else NA,
        pval = if (inf_obj_res$status == "Success") inf_obj_res$pval else NA,
        bootstrap_ci_lower = if (inf_obj_res$status == "Success") inf_obj_res$bootstrap_ci_lower else NA_real_,
        bootstrap_ci_upper = if (inf_obj_res$status == "Success") inf_obj_res$bootstrap_ci_upper else NA_real_,
        bootstrap_ci_error = if (inf_obj_res$status == "Success") inf_obj_res$bootstrap_ci_error else NA_character_,
        bootstrap_pval = if (inf_obj_res$status == "Success") inf_obj_res$bootstrap_pval else NA_real_,
        bootstrap_pval_error = if (inf_obj_res$status == "Success") inf_obj_res$bootstrap_pval_error else NA_character_,
        rand_pval = if (inf_obj_res$status == "Success") inf_obj_res$rand_pval else NA_real_,
        rand_pval_error = if (inf_obj_res$status == "Success") inf_obj_res$rand_pval_error else NA_character_,
        rand_ci_lower = if (inf_obj_res$status == "Success") inf_obj_res$rand_ci_lower else NA_real_,
        rand_ci_upper = if (inf_obj_res$status == "Success") inf_obj_res$rand_ci_upper else NA_real_,
        rand_ci_error = if (inf_obj_res$status == "Success") inf_obj_res$rand_ci_error else NA_character_,
        mle_ci_lower = if (inf_obj_res$status == "Success") inf_obj_res$mle_ci_lower else NA_real_,
        mle_ci_upper = if (inf_obj_res$status == "Success") inf_obj_res$mle_ci_upper else NA_real_,
        mle_ci_error = if (inf_obj_res$status == "Success") inf_obj_res$mle_ci_error else NA_character_,
        mle_pval = if (inf_obj_res$status == "Success") inf_obj_res$mle_pval else NA_real_,
        mle_pval_error = if (inf_obj_res$status == "Success") inf_obj_res$mle_pval_error else NA_character_
      )
    }
  }
}

final_results <- rbindlist(results)
print(final_results[status == "Success", .(
  response_type,
  design,
  inference,
  time_ms,
  estimate,
  pval,
  bootstrap_ci_lower,
  bootstrap_ci_upper,
  bootstrap_pval,
  rand_ci_lower,
  rand_ci_upper,
  rand_pval,
  mle_ci_lower,
  mle_ci_upper,
  mle_pval
)])

# Summary of errors
errors <- final_results[status != "Success", .(count = .N), by = .(status)]
if (nrow(errors) > 0) {
  vbcat("\nInference compatibility/error summary:\n")
  print(errors)
}

fwrite(final_results, "benchmark_results.csv")
vbcat("\nBenchmark complete. Results saved to benchmark_results.csv\n")
