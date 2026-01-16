if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("benchmark_inference_loops.R requires the data.table package")
}

library(data.table)

safe_run_timed <- function(expr) {
  start <- proc.time()
  result <- tryCatch({
    list(value = expr(), error = NA_character_)
  }, error = function(e) {
    list(value = NULL, error = conditionMessage(e))
  })
  result$elapsed <- as.numeric((proc.time() - start)[["elapsed"]])
  result
}

designs <- c("SeqDesignCRD", "SeqDesignEfron", "SeqDesignAtkinson", "SeqDesignKK14", "SeqDesignKK21", "SeqDesignKK21stepwise", "SeqDesigniBCRD")

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

generate_responses <- function(n) {
  list(
    continuous = list(y = rnorm(n), dead = rep(1, n)),
    incidence = list(y = rbinom(n, 1, 0.5), dead = rep(1, n)),
    count = list(y = rpois(n, 5), dead = rep(1, n)),
    proportion = list(y = runif(n, 0.01, 0.99), dead = rep(1, n)),
    survival = list(y = rexp(n, 0.1), dead = rbinom(n, 1, 0.8))
  )
}

run_comprehensive_inference_suite <- function(
  version_label,
  iteration = NA_integer_,
  n = 100,
  p = 5,
  seed = 123,
  verbose = TRUE
) {
  set.seed(seed)
  responses <- generate_responses(n)
  X_mat <- matrix(rnorm(n * p), nrow = n)
  colnames(X_mat) <- paste0("x", 1:p)
  X_dt <- as.data.table(X_mat)

  results <- list()
  row_idx <- 0

  for (rt in names(inf_by_res_type)) {
    resp_info <- responses[[rt]]
    y <- resp_info$y
    dead <- resp_info$dead

    for (des_name in designs) {
      des_class <- get(des_name, inherits = TRUE)
      des_init_res <- safe_run_timed(function() {
        des_obj <- des_class$new(n = n, response_type = rt, verbose = verbose)
        for (i in seq_len(n)) {
          des_obj$add_subject_to_experiment_and_assign(X_dt[i, ])
        }
        des_obj$add_all_subject_responses(y, deads = dead)
        des_obj
      })

      row_idx <- row_idx + 1
      results[[row_idx]] <- data.table(
        version_label = version_label,
        iteration = iteration,
        response_type = rt,
        design = des_name,
        inference = NA_character_,
        metric = "design_initialization",
        elapsed = des_init_res$elapsed,
        error = des_init_res$error,
        status = ifelse(is.na(des_init_res$error), "success", "error")
      )

      if (!is.na(des_init_res$error)) {
        next
      }
      des_obj <- des_init_res$value

      inf_methods <- inf_by_res_type[[rt]]
      for (inf_name in inf_methods) {
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
          row_idx <- row_idx + 1
          results[[row_idx]] <- data.table(
            version_label = version_label,
            iteration = iteration,
            response_type = rt,
            design = des_name,
            inference = inf_name,
            metric = "skip",
            elapsed = NA_real_,
            error = skip_reason,
            status = skip_reason
          )
          next
        }

        inf_class <- get(inf_name, inherits = TRUE)
        init_res <- safe_run_timed(function() inf_class$new(des_obj, verbose = verbose))
        row_idx <- row_idx + 1
        results[[row_idx]] <- data.table(
          version_label = version_label,
          iteration = iteration,
          response_type = rt,
          design = des_name,
          inference = inf_name,
          metric = "initialize",
          elapsed = init_res$elapsed,
          error = init_res$error,
          status = ifelse(is.na(init_res$error), "success", "error")
        )

        if (!is.na(init_res$error)) {
          next
        }
        inf_obj <- init_res$value

        metrics <- list(
          list(name = "compute_treatment_estimate", fn = function() inf_obj$compute_treatment_estimate()),
          list(name = "compute_bootstrap_confidence_interval", fn = function() inf_obj$compute_bootstrap_confidence_interval(B = 100)),
          list(name = "compute_bootstrap_two_sided_pval", fn = function() inf_obj$compute_bootstrap_two_sided_pval(B = 100)),
          list(name = "compute_two_sided_pval_for_treatment_effect_rand", fn = function() inf_obj$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = 100)),
          list(name = "compute_confidence_interval_rand", fn = function() inf_obj$compute_confidence_interval_rand(nsim_exact_test = 100))
        )

        if ("compute_mle_confidence_interval" %in% names(inf_obj)) {
          metrics[[length(metrics) + 1]] <- list(
            name = "compute_mle_confidence_interval",
            fn = function() inf_obj$compute_mle_confidence_interval()
          )
        }
        if ("compute_mle_two_sided_pval_for_treatment_effect" %in% names(inf_obj)) {
          metrics[[length(metrics) + 1]] <- list(
            name = "compute_mle_two_sided_pval_for_treatment_effect",
            fn = function() inf_obj$compute_mle_two_sided_pval_for_treatment_effect()
          )
        }

        for (metric in metrics) {
          metric_res <- safe_run_timed(metric$fn)
          row_idx <- row_idx + 1
          results[[row_idx]] <- data.table(
            version_label = version_label,
            iteration = iteration,
            response_type = rt,
            design = des_name,
            inference = inf_name,
            metric = metric$name,
            elapsed = metric_res$elapsed,
            error = metric_res$error,
            status = ifelse(is.na(metric_res$error), "success", "error")
          )
        }
      }
    }
  }

  data.table::rbindlist(results, fill = TRUE)
}
