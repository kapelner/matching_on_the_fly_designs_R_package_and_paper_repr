#' Asymptotic Inference
#'
#' Abstract class for asymptotic inference.
#'
#' @keywords internal
InferenceAsymp = R6::R6Class("InferenceAsymp",
	lock_objects = FALSE,
	inherit = InferenceBoot,
	public = list(
		#' @description
		#' Computes an asymptotic confidence interval based on Wald-type tests.
		#'
		#' @param alpha					Significance level 1 - \code{alpha}. Default 0.05.
		#'
		#' @return 	A Wald-type confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			
			est = self$compute_treatment_estimate()
			se = private$get_standard_error()
			df = private$get_degrees_of_freedom()
			
			if (!is.finite(est) || !is.finite(se) || se <= 0) return(c(NA_real_, NA_real_))
			
			critical_val = if (is.finite(df)) stats::qt(1 - alpha / 2, df = df) else stats::qnorm(1 - alpha / 2)
			
			ci = c(est - critical_val * se, est + critical_val * se)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		#' @description
		#' Computes an asymptotic two-sided p-value for the treatment effect.
		#'
		#' @param delta					Null treatment effect to test against. Default 0.
		#'
		#' @return 	The asymptotic p-value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			
			est = self$compute_treatment_estimate()
			se = private$get_standard_error()
			df = private$get_degrees_of_freedom()
			
			if (!is.finite(est) || !is.finite(se) || se <= 0) return(NA_real_)
			
			t_stat = (est - delta) / se
			
			if (is.finite(df)) {
				2 * stats::pt(-abs(t_stat), df = df)
			} else {
				2 * stats::pnorm(-abs(t_stat))
			}
		},

		#' @description
		#' Abstract method to compute the treatment estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		#' @return 	A scalar treatment estimate.
		compute_treatment_estimate = function(estimate_only = FALSE){
			stop("Must be implemented by concrete class.")
		},

		#' @description
		#' Returns the model object from the last call that produced the treatment
		#' estimate and SE. Calls \code{compute_treatment_estimate()} first if needed.
		#'
		#' @return The cached model object (type depends on the concrete class).
		get_mod = function(){
			if (is.null(private$cached_mod)) self$compute_treatment_estimate()
			private$cached_mod
		},

		#' @description
		#' Prints a summary of the model from the last call that produced the
		#' treatment estimate and SE.
		get_summary = function(){
			mod = self$get_mod()
			if (is.null(mod)) {
				cat("No model available (call compute_treatment_estimate() first).\n")
				return(invisible(NULL))
			}
			if (identical(class(mod), "list")) {
				if (!is.null(private$cached_values$summary_table)) {
					print(private$cached_values$summary_table)
				} else {
					print(mod)
				}
			} else {
				print(summary(mod))
			}
			invisible(NULL)
		}
	),

	private = list(
		cached_mod = NULL,
		get_standard_error = function() stop("Must be implemented by concrete class or shared helper."),
		get_degrees_of_freedom = function() NA_real_,

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
		
		# Shared helpers for z/t tests
		compute_z_or_t_ci_from_s_and_df = function(alpha){
			beta_hat_T = private$cached_values$beta_hat_T
			s_beta_hat_T = private$cached_values$s_beta_hat_T
			is_z = private$cached_values$is_z
			df = private$cached_values$df
			
			if (length(beta_hat_T) != 1L || length(s_beta_hat_T) != 1L) return(c(NA_real_, NA_real_))
			if (!is.finite(beta_hat_T) || !is.finite(s_beta_hat_T) || s_beta_hat_T <= 0) return(c(NA_real_, NA_real_))
			
			mult = if (isTRUE(is_z) || !is.finite(df)) stats::qnorm(1 - alpha / 2) else stats::qt(1 - alpha / 2, df = df)
			ci = c(beta_hat_T - mult * s_beta_hat_T, beta_hat_T + mult * s_beta_hat_T)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		compute_z_or_t_two_sided_pval_from_s_and_df = function(delta){
			beta_hat_T = private$cached_values$beta_hat_T
			s_beta_hat_T = private$cached_values$s_beta_hat_T
			is_z = private$cached_values$is_z
			df = private$cached_values$df
			
			if (length(beta_hat_T) != 1L || length(s_beta_hat_T) != 1L) return(NA_real_)
			if (!is.finite(beta_hat_T) || !is.finite(s_beta_hat_T) || s_beta_hat_T <= 0) return(NA_real_)
			
			val = (beta_hat_T - delta) / s_beta_hat_T
			if (isTRUE(is_z) || !is.finite(df)) 2 * stats::pnorm(-abs(val)) else 2 * stats::pt(-abs(val), df = df)
		}
	)
)
