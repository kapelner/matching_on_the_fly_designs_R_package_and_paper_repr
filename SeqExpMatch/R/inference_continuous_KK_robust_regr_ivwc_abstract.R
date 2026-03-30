#' Abstract class for Robust-Regression IVWC Compound Inference for KK Designs
#'
#' Fits a variance-weighted compound estimator for KK matching-on-the-fly designs
#' with continuous responses using robust linear regression (`MASS::rlm`) for the
#' matched-pair and reservoir components separately.
#'
#' @keywords internal
InferenceAbstractKKRobustRegrIVWC = R6::R6Class("InferenceAbstractKKRobustRegrIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThroughCompound,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param method			Robust-regression fitting method for `MASS::rlm`; one of `"M"` or `"MM"`.
		#' @param num_cores			Number of CPU cores for parallel processing.
		#' @param verbose			Whether to print progress messages.
		#' @param make_fork_cluster Whether to use a fork cluster for parallelization.
		initialize = function(des_obj, method = "MM", num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "continuous")
			assertChoice(method, c("M", "MM"))
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
			private$rlm_method = method
		},

		#' @description
		#' Returns the estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the approximate confidence interval.
		#' @param alpha The confidence level in the computed confidence interval is 1 -
		#'   \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the approximate p-value.
		#' @param delta The null difference to test against. For any treatment effect at all this
		#'   is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		#' @description
		#' Duplicate
		#' @param verbose A flag indicating whether messages should be displayed.
		duplicate = function(verbose = FALSE){
			i = super$duplicate(verbose = verbose)
			i
		}
	),

	private = list(
		rlm_method = NULL,
		rlm_force_M = FALSE,

		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		reduce_design_matrix_once = function(X, j_treat, cache_key){
			cached = private$cached_values[[cache_key]]
			if (!is.null(cached)) return(cached)

			qr_X = qr(X)
			if (qr_X$rank < ncol(X)){
				keep = qr_X$pivot[seq_len(qr_X$rank)]
				if (!(j_treat %in% keep)) keep[qr_X$rank] = j_treat
				keep = sort(unique(keep))
				X = X[, keep, drop = FALSE]
				j_treat = which(keep == j_treat)
			}

			cached = list(X = X, j_treat = j_treat)
			private$cached_values[[cache_key]] = cached
			cached
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			if (m > 0){
				private$robust_for_matched_pairs()
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m  = private$cached_values$ssq_beta_T_matched
			m_ok   = !is.null(beta_m) && is.finite(beta_m) && !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0

			if (nRT > 0 && nRC > 0){
				private$robust_for_reservoir()
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r  = private$cached_values$ssq_beta_T_reservoir
			r_ok   = !is.null(beta_r) && is.finite(beta_r) && !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0

			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T   = beta_m
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T   = beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		fit_rlm_with_treatment = function(X, y, j_treat){
			if (nrow(X) <= ncol(X)) return(NULL)

			run_rlm = function(method){
				tryCatch({
					if (identical(method, "M")) {
						MASS::rlm(x = X, y = y, method = "M", psi = MASS::psi.huber)
					} else {
						MASS::rlm(x = X, y = y, method = method)
					}
				}, error = function(e) e)
			}

			method_to_try = if (isTRUE(private$rlm_force_M)) "M" else private$rlm_method
			mod = run_rlm(method_to_try)
			if (inherits(mod, "error") && identical(method_to_try, "MM")){
				msg = if (length(mod$message) == 0L) "" else mod$message
				if (grepl("'lqs' failed", msg, fixed = TRUE) || grepl("singular", msg, ignore.case = TRUE)) {
					private$rlm_force_M = TRUE
					mod = run_rlm("M")
				}
			}
			if (inherits(mod, "error")) return(NULL)
			if (is.null(mod)) return(NULL)

			coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			if ((is.null(coef_table) || nrow(coef_table) < j_treat) && identical(method_to_try, "MM")){
				private$rlm_force_M = TRUE
				mod = run_rlm("M")
				if (inherits(mod, "error") || is.null(mod)) return(NULL)
				coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			}
			if (is.null(coef_table) || nrow(coef_table) < j_treat) return(NULL)

			beta = as.numeric(coef_table[j_treat, "Value"])
			se   = as.numeric(coef_table[j_treat, "Std. Error"])
			if ((!is.finite(beta) || !is.finite(se) || se <= 0) && identical(method_to_try, "MM")){
				private$rlm_force_M = TRUE
				mod = run_rlm("M")
				if (inherits(mod, "error") || is.null(mod)) return(NULL)
				coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
				if (is.null(coef_table) || nrow(coef_table) < j_treat) return(NULL)
				beta = as.numeric(coef_table[j_treat, "Value"])
				se   = as.numeric(coef_table[j_treat, "Std. Error"])
			}
			if (!is.finite(beta) || !is.finite(se) || se <= 0) return(NULL)

			list(beta = beta, ssq = se^2)
		},

		robust_for_matched_pairs = function(){
			yd = private$cached_values$KKstats$y_matched_diffs
			m  = length(yd)
			if (private$include_covariates()){
				Xd = as.matrix(private$cached_values$KKstats$X_matched_diffs)
				X = if (ncol(Xd) > 0L) cbind(1, Xd) else matrix(1, nrow = m, ncol = 1L)
				reduced = private$reduce_design_matrix_once(
					X,
					1L,
					cache_key = "kk_robust_ivwc_matched_reduced_design"
				)
				X = reduced$X
			} else {
				X = matrix(1, nrow = m, ncol = 1L)
			}

			fit = private$fit_rlm_with_treatment(X, yd, 1L)
			if (is.null(fit)) {
				private$cached_values$beta_T_matched     = if (m >= 1) mean(yd) else NA_real_
				private$cached_values$ssq_beta_T_matched = if (m >= 2) var(yd) / m else NA_real_
			} else {
				private$cached_values$beta_T_matched     = fit$beta
				private$cached_values$ssq_beta_T_matched = fit$ssq
			}
		},

		robust_for_reservoir = function(){
			y_r = private$cached_values$KKstats$y_reservoir
			w_r = private$cached_values$KKstats$w_reservoir
			X_r = as.matrix(private$cached_values$KKstats$X_reservoir)
			j_treat = 2L

			if (private$include_covariates()){
				X_full = cbind(1, w_r, X_r)
				reduced = private$reduce_design_matrix_once(
					X_full,
					j_treat,
					cache_key = "kk_robust_ivwc_reservoir_reduced_design"
				)
				X_full = reduced$X
				j_treat = reduced$j_treat
			} else {
				X_full = cbind(1, w_r)
			}

			fit = private$fit_rlm_with_treatment(X_full, y_r, j_treat)
			if (is.null(fit)) {
				private$cached_values$beta_T_reservoir     = NA_real_
				private$cached_values$ssq_beta_T_reservoir = NA_real_
			} else {
				private$cached_values$beta_T_reservoir     = fit$beta
				private$cached_values$ssq_beta_T_reservoir = fit$ssq
			}
		}
	)
)
