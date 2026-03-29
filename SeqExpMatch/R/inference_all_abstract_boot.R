#' Bootstrap-based Inference
#'
#' Abstract class for bootstrap-based inference.
#'
#' @keywords internal
InferenceBoot = R6::R6Class("InferenceBoot",
	lock_objects = FALSE,
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

			# Duplicate objects for thread safety
			inf_template = self$duplicate()
			des_template = private$des_obj$duplicate()
			is_KK_local = private$is_KK

			# Determine cores — warm-up guard: run one iteration serially to estimate per-iteration
			# cost, then only parallelize if computation outweighs overhead per worker.
			# For a persistent fork cluster the per-call overhead is ~10ms (socket round-trip);
			# for mclapply (per-call fork) it is ~500ms.
			actual_cores = private$num_cores
			if (actual_cores > 1L) {
				t_boot_warmup = system.time({
					w_des = des_template$duplicate()
					w_inf = inf_template$duplicate(); w_inf$.__enclos_env__$private$num_cores = 1L
					w_des$resample_design()
					w_inf$.__enclos_env__$private$w = w_des$.__enclos_env__$private$w
					w_inf$.__enclos_env__$private$y = w_des$.__enclos_env__$private$y
					if (is_KK_local && !is.null(w_inf$.__enclos_env__$private$compute_basic_match_data)) {
						w_inf$.__enclos_env__$private$m = w_des$.__enclos_env__$private$m
						w_inf$.__enclos_env__$private$compute_basic_match_data()
					}
					tryCatch(w_inf$compute_treatment_estimate(estimate_only = TRUE), error = function(e) NA_real_)
				})[[3]]
				fork_overhead_estimate = if (isTRUE(private$make_fork_cluster)) 0.01 else 0.5
				if (!(t_boot_warmup * B > fork_overhead_estimate * actual_cores))
					actual_cores = 1L
			}

			# Fork-cluster path: use parLapply with closure capture.
			# (clusterExport + environment(fn)=globalenv() is broken: parLapply serializes
			# a snapshot of master's globalenv, not the worker's live globalenv.)
			if (isTRUE(private$make_fork_cluster) && actual_cores > 1L) {
				cl = private$get_or_create_fork_cluster()
				boot_distr = unlist(parallel::parLapply(cl, 1:B, function(idx) {
					try(set_package_threads(1L), silent = TRUE)
					worker_des = des_template$duplicate()
					worker_inf = inf_template$duplicate()
					worker_inf$.__enclos_env__$private$num_cores = 1L
					worker_des$resample_design()
					worker_inf$.__enclos_env__$private$w = worker_des$.__enclos_env__$private$w
					worker_inf$.__enclos_env__$private$y = worker_des$.__enclos_env__$private$y
					if (is_KK_local && !is.null(worker_inf$.__enclos_env__$private$compute_basic_match_data)) {
						worker_inf$.__enclos_env__$private$m = worker_des$.__enclos_env__$private$m
						worker_inf$.__enclos_env__$private$compute_basic_match_data()
					}
					tryCatch(worker_inf$compute_treatment_estimate(estimate_only = TRUE), error = function(e) NA_real_)
				}))
			} else {
				boot_distr = unlist(private$par_lapply(1:B, function(idx) {
					set_package_threads(1L)
					worker_des = des_template$duplicate()
					worker_inf = inf_template$duplicate()
					worker_inf$.__enclos_env__$private$num_cores = 1L
					worker_des$resample_design()
					worker_inf$.__enclos_env__$private$w = worker_des$.__enclos_env__$private$w
					worker_inf$.__enclos_env__$private$y = worker_des$.__enclos_env__$private$y
					if (is_KK_local && !is.null(worker_inf$.__enclos_env__$private$compute_basic_match_data)) {
						worker_inf$.__enclos_env__$private$m = worker_des$.__enclos_env__$private$m
						worker_inf$.__enclos_env__$private$compute_basic_match_data()
					}
					tryCatch(worker_inf$compute_treatment_estimate(estimate_only = TRUE), error = function(e) NA_real_)
				}, n_cores = actual_cores, show_progress = show_progress))
			}

			if (!is.numeric(boot_distr)) boot_distr = as.numeric(boot_distr)
			
			if (is.null(private$cached_values$boot_distr_cache)) private$cached_values$boot_distr_cache = list()
			private$cached_values$boot_distr_cache[[cache_key]] = boot_distr
			boot_distr
		},

		#' @description
		#' Computes a bootstrap-based two-sided p-value for the treatment effect.
		#'
		#' @param delta					Null hypothesis value. Default 0.
		#' @param B						Number of bootstrap samples. Default 501.
		#' @param na.rm					Remove non-finite bootstrap replicates. Default FALSE.
		#'
		#' @return 	A bootstrap two-sided p-value.
		compute_bootstrap_two_sided_pval = function(delta = 0, B = 501, na.rm = FALSE){
			assertNumeric(delta, len = 1)
			assertCount(B, positive = TRUE)
			assertFlag(na.rm)
			boot_distr = self$approximate_bootstrap_distribution_beta_hat_T(B)
			if (isTRUE(na.rm)) boot_distr = boot_distr[is.finite(boot_distr)]
			else if (any(!is.finite(boot_distr))) return(NA_real_)
			if (length(boot_distr) == 0L) return(NA_real_)
			est = as.numeric(self$compute_treatment_estimate())
			if (length(est) == 0L || !is.finite(est[1])) return(NA_real_)
			est = est[1]
			# Shift bootstrap distribution to be centered at delta (null hypothesis)
			boot_null = boot_distr - mean(boot_distr) + delta
			n_bs = length(boot_null)
			min(1, max(2 / n_bs, 2 * min(
				sum(boot_null >= est) / n_bs,
				sum(boot_null <= est) / n_bs
			)))
		},

		#' @description
		#' Computes a bootstrap-based confidence interval.
		#'
		#' @param alpha					The confidence level 1 - \code{alpha}. Default 0.05.
		#' @param B						Number of bootstrap samples. Default 501.
		#' @param type					Type of bootstrap CI: "percentile" (default) or "basic".
		#' @param na.rm                                   Remove non-finite bootstrap replicates.
		#'   Default TRUE. Non-finite replicates are always removed internally.
		#' @param show_progress			Show progress bar.
		#'
		#' @return 	A bootstrap confidence interval.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = "percentile", na.rm = TRUE, show_progress = TRUE){
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
