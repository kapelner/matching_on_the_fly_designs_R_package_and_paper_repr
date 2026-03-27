# Bootstrap-based Inference
#
# @description
# Abstract class for bootstrap-based inference.
#
# @keywords internal
InferenceBoot = R6::R6Class("InferenceBoot",
	inherit = InferenceRandCI,
	public = list(
		#' @description
		#' Creates the bootstrap distribution of the estimate for the treatment effect
		#'
		#' @param B						Number of bootstrap samples. The default is 501.
		#' @param show_progress			A flag indicating whether a progress bar should be displayed.
		#'
		#' @return 	A vector of length \code{B} containing the bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE){
			private$assert_design_supports_resampling("Bootstrap inference")
			assertCount(B, positive = TRUE)

			# Check cache
			cache_key = as.character(B)
			if (!is.null(private$cached_values$boot_distr_cache[[cache_key]])) {
				return(private$cached_values$boot_distr_cache[[cache_key]])
			}

			if (private$verbose) cat("Computing bootstrap distribution...\n")

			mclapply_fn = if (isTRUE(show_progress) && requireNamespace("pbmcapply", quietly = TRUE)) pbmcapply::pbmclapply else parallel::mclapply

			# Duplicate objects for thread safety
			inf_template = self$duplicate()
			des_template = private$des_obj$duplicate()

			# Determine cores
			actual_cores = private$num_cores
			
			boot_distr = unlist(mclapply_fn(1:B, function(idx) {
				set_package_threads(1L)
				worker_des = des_template$duplicate()
				worker_inf = inf_template$duplicate()
				worker_inf$.__enclos_env__$private$num_cores = 1L
				
				# Resample design
				worker_des$resample_design()
				
				# Update inference object with resampled design
				worker_inf$.__enclos_env__$private$w = worker_des$.__enclos_env__$private$w
				worker_inf$.__enclos_env__$private$y = worker_des$.__enclos_env__$private$y
				if (private$is_KK) {
					worker_inf$.__enclos_env__$private$m = worker_des$.__enclos_env__$private$m
					worker_inf$.__enclos_env__$private$compute_basic_match_data()
				}
				
				# Compute estimate
				tryCatch(worker_inf$compute_treatment_estimate(), error = function(e) NA_real_)
			}, mc.cores = actual_cores))

			if (!is.numeric(boot_distr)) boot_distr = as.numeric(boot_distr)
			
			if (is.null(private$cached_values$boot_distr_cache)) private$cached_values$boot_distr_cache = list()
			private$cached_values$boot_distr_cache[[cache_key]] = boot_distr
			boot_distr
		},

		#' @description
		#' Computes a bootstrap-based confidence interval.
		#'
		#' @param alpha					The confidence level 1 - \code{alpha}. Default 0.05.
		#' @param B						Number of bootstrap samples. Default 501.
		#' @param type					Type of bootstrap CI: "percentile" (default) or "basic".
		#' @param show_progress			Show progress bar.
		#'
		#' @return 	A bootstrap confidence interval.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = "percentile", show_progress = TRUE){
			private$assert_design_supports_resampling("Bootstrap inference")
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertChoice(type, c("percentile", "basic"))

			boot_distr = self$approximate_bootstrap_distribution_beta_hat_T(B, show_progress)
			boot_distr = boot_distr[is.finite(boot_distr)]
			if (length(boot_distr) < B / 2) stop("Too many bootstrap samples failed.")

			est = self$compute_treatment_estimate()
			
			ci = if (type == "percentile") {
				stats::quantile(boot_distr, probs = c(alpha / 2, 1 - alpha / 2))
			} else {
				2 * est - stats::quantile(boot_distr, probs = c(1 - alpha / 2, alpha / 2))
			}
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	),
	private = list(
		# Cache for bootstrap distributions
		boot_distr_cache = list()
	)
)
