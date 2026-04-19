#' Robust Regression Inference for Continuous Responses
#'
#' Fits a robust linear regression via \code{MASS::rlm} for continuous responses
#' using the treatment indicator and, optionally, all recorded covariates as
#' predictors. This provides a Huber/MM-style robustness upgrade over ordinary
#' least squares when outcomes are heavy-tailed or outlier-prone. Inference is
#' based on the coefficient table returned by \code{summary.rlm()}.
#'
#' @details
#' The \code{method} argument is passed to \code{MASS::rlm} and may be either
#' \code{"M"} or \code{"MM"}. For \code{"M"}, the fit uses Huber's psi
#' function. Approximate confidence intervals and p-values use the reported
#' robust standard error with residual degrees of freedom \eqn{n - p}.
#'
#' @export
InferenceContinRobustRegr = R6::R6Class("InferenceContinRobustRegr",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize a robust-regression inference object for a completed design
		#' with a continuous response.
		#' @param des_obj A completed \code{Design} object with a continuous response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates and \code{FALSE} otherwise.
		#' @param method Robust-regression fitting method for \code{MASS::rlm}; one
		#'   of \code{"M"} or \code{"MM"}. The default is \code{"MM"}.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use the optimized Rcpp
		#'   implementation. If \code{FALSE}, use \code{MASS::rlm}.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, method = "MM", use_rcpp = TRUE, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "continuous")
				assertChoice(method, c("M", "MM"))
				assertFlag(include_covariates, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			super$initialize(des_obj, verbose = verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			
			if (is.null(include_covariates)) {
				include_covariates = des_obj$has_covariates()
			}
			private$include_covariates = include_covariates
			private$rlm_method = method
			private$use_rcpp = use_rcpp
		},

		#' @description
		#' Computes the robust-regression estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes an approximate confidence interval for the treatment effect.
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes an approximate two-sided p-value for the treatment effect.
		#' @param delta The null difference to test against. Default is zero.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		rlm_method = NULL,
		use_rcpp = TRUE,
		include_covariates = NULL,
		fit_warm_coefficients = NULL,
		fit_warm_keep = NULL,
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_treatment_estimate(estimate_only = estimate_only))
			}

			# Use the same design matrix structure as the original fit
			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()

			if (length(Xmm_cols) == 0L){
				# Univariate case
				X_fit = cbind(1, treatment = private$w)
			} else {
				# Multivariate case
				X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				X_fit = cbind(1, treatment = private$w, X_cov)
			}

			if (private$use_rcpp) {
				fit = private$fit_rlm_model(X_fit, estimate_only = estimate_only, warm_start = TRUE)
				if (is.null(fit)) return(NA_real_)
				return(as.numeric(fit$coefficients[2]))
			} else {
				fit = tryCatch(
					suppressWarnings(MASS::rlm(x = X_fit, y = as.numeric(private$y), method = private$rlm_method, scale.est = "mad", maxit = 20)),
					error = function(e) NULL
				)
				if (is.null(fit) || !is.finite(stats::coef(fit)[2])){
					return(NA_real_)
				}
				return(as.numeric(stats::coef(fit)[2]))
			}
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		create_bootstrap_worker_state = function(){
			private$create_design_backed_bootstrap_worker_state()
		},

		load_bootstrap_sample_into_worker = function(worker_state, indices){
			private$load_bootstrap_sample_into_design_backed_worker(worker_state, indices)
		},

		compute_bootstrap_worker_estimate = function(worker_state){
			private$compute_bootstrap_worker_estimate_via_compute_treatment_estimate(worker_state)
		},

		build_design_matrix = function(){
			if (private$include_covariates) {
				private$create_design_matrix()
			} else {
				X = cbind(1, private$w)
				colnames(X) = c("(Intercept)", "treatment")
				X
			}
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},

		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df
		},

		set_failed_fit_cache = function(){
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$is_z = FALSE
			private$cached_values$df = NA_real_
		},

		get_ci_fit_controls = function(){
			ctrl = private$randomization_mc_control
			list(
				warm_start = !is.null(ctrl) && isTRUE(ctrl$fit_warm_start_enable),
				reuse_factorizations = !is.null(ctrl) && isTRUE(ctrl$fit_reuse_factorizations)
			)
		},

		fit_rlm_model = function(X_fit, estimate_only = FALSE, warm_start = FALSE){
			start_beta = if (warm_start && !is.null(private$fit_warm_coefficients) && 
							 length(private$fit_warm_coefficients) == ncol(X_fit)) private$fit_warm_coefficients else NULL
			
			if (private$use_rcpp) {
				tryCatch(
					fast_robust_regression_cpp(X = X_fit, y = as.numeric(private$y), start_beta = start_beta, method = private$rlm_method, j = 2L),
					error = function(e) NULL
				)
			} else {
				tryCatch(
					suppressWarnings(MASS::rlm(x = X_fit, y = as.numeric(private$y), method = private$rlm_method, init = start_beta)),
					error = function(e) NULL
				)
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			fit_controls = private$get_ci_fit_controls()
			
			if (is.null(private$best_Xmm_colnames)) {
				X_full = private$build_design_matrix()
				attempt = private$fit_with_hardened_qr_column_dropping(
					X_full = X_full,
					required_cols = match("treatment", colnames(X_full)),
					fit_fun = function(X_fit, j_treat){
						mod = private$fit_rlm_model(X_fit, estimate_only = estimate_only, warm_start = fit_controls$warm_start)
						if (is.null(mod)) return(NULL)
						mod$j_treat = j_treat
						mod
					},
					fit_ok = function(mod, X_fit, keep){
						j_treat = mod$j_treat
						beta_hat = as.numeric(if (private$use_rcpp) mod$coefficients else stats::coef(mod))
						if (is.null(mod) || length(beta_hat) < j_treat || !is.finite(beta_hat[j_treat])) return(FALSE)
						if (estimate_only) return(TRUE)
						
						if (private$use_rcpp) {
							se_w = as.numeric(sqrt(mod$ssq_b_j))
						} else {
							st = tryCatch(summary(mod), error = function(e) NULL)
							if (is.null(st)) return(FALSE)
							se_w = as.numeric(st$coefficients[j_treat, "Std. Error"])
						}
						is.finite(se_w) && se_w > 0
					}
				)

				if (is.null(attempt$fit)){
					private$set_failed_fit_cache()
					return(invisible(NULL))
				}
				
				fit = attempt$fit
				X_fit = attempt$X_fit
				j_treat = fit$j_treat
				private$best_Xmm_colnames = setdiff(colnames(X_fit), c("(Intercept)", "treatment"))
				private$fit_warm_coefficients = as.numeric(if (private$use_rcpp) fit$coefficients else stats::coef(fit))
			} else {
				# Reuse structure
				X_data = private$get_X()
				Xmm_cols = private$best_Xmm_colnames
				if (length(Xmm_cols) == 0L){
					X_fit = cbind(1, treatment = private$w)
				} else {
					X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
					X_fit = cbind(1, treatment = private$w, X_cov)
				}
				j_treat = 2L
				fit = private$fit_rlm_model(X_fit, estimate_only = estimate_only, warm_start = fit_controls$warm_start)
				if (is.null(fit)) {
					private$set_failed_fit_cache()
					return(invisible(NULL))
				}
			}

			beta_hat = as.numeric(if (private$use_rcpp) fit$coefficients else stats::coef(fit))
			private$cached_values$beta_hat_T = beta_hat[j_treat]
			if (estimate_only) return(invisible(NULL))

			if (private$use_rcpp) {
				private$cached_values$s_beta_hat_T = as.numeric(sqrt(fit$ssq_b_j))
			} else {
				st = summary(fit)
				private$cached_values$s_beta_hat_T = as.numeric(st$coefficients[j_treat, "Std. Error"])
			}
			private$cached_values$df = nrow(X_fit) - ncol(X_fit)
			private$cached_values$is_z = FALSE # rlm uses t-distribution by default
		}
	)
)
