#' Hurdle Negative Binomial Inference for Count Responses
#'
#' Internal base class for non-KK hurdle negative binomial regression models. The
#' hurdle indicator model is fit on all subjects, and the count component is fit on
#' the positive-count subjects using a zero-truncated negative binomial model
#' optimized in C++ via L-BFGS. The reported treatment effect is the treatment
#' coefficient from the conditional count component, on the log-rate scale.
#'
#' @keywords internal
#' @noRd
InferenceCountHurdleNegBinAbstract = R6::R6Class("InferenceCountHurdleNegBinAbstract",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param verbose A flag indicating whether messages should be displayed.
		initialize = function(des_obj, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "count")
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute asymp confidence interval
		#' @param alpha Description for alpha
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute asymp two sided pval for treatment effect
		#' @param delta Description for delta
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		hurdle_description = function() stop(class(self)[1], " must implement hurdle_description()."),

		predictors_df = function(){
			data.frame(w = private$w)
		},

		try_hurdle_negbin_fit = function(X_f, j_t, estimate_only = FALSE){
			mod = tryCatch(
				if (estimate_only) {
					fast_hurdle_negbin_cpp(X_f, private$y)
				} else {
					fast_hurdle_negbin_with_var_cpp(X_f, private$y, j = j_t)
				},
				error = function(e) NULL
			)
			if (is.null(mod)) return(NULL)
			b = as.numeric(mod$b)
			ssq = if (estimate_only) NA_real_ else as.numeric(mod$ssq_b_j)
			if (length(b) != ncol(X_f) || any(!is.finite(b))) return(NULL)
			if (!estimate_only && (!is.finite(ssq) || ssq < 0)) return(NULL)
			list(mod = mod, b = b, ssq = ssq, j = j_t, X = X_f)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			X_full = cbind(1, private$w, as.matrix(private$predictors_df()[, setdiff(colnames(private$predictors_df()), "w"), drop = FALSE]))
			colnames(X_full)[1:2] = c("(Intercept)", "treatment")
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$cache_nonestimable_estimate("hurdle_negbin_design_unusable")
				return(invisible(NULL))
			}

			fit = private$try_hurdle_negbin_fit(X_fit, j_treat, estimate_only = estimate_only)
			if (private$harden && is.null(fit) && ncol(X_fit) > 2L){
				X_treat_only = X_fit[, 1:2, drop = FALSE]
				fit = private$try_hurdle_negbin_fit(X_treat_only, 2L, estimate_only = estimate_only)
				reduced = list(X = X_treat_only, keep = 1:2, j_treat = 2L)
			}
			if (is.null(fit)){
				private$cache_nonestimable_estimate("hurdle_negbin_fit_unavailable")
				return(invisible(NULL))
			}

			b_full = rep(NA_real_, ncol(X_full))
			b_full[reduced$keep] = fit$b
			names(b_full) = colnames(X_full)

			private$cached_values$beta_hat_T = fit$b[fit$j]
			if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(fit$ssq)
			private$cached_values$is_z = TRUE
			private$cached_values$full_coefficients = b_full
			private$cached_values$theta_hat = as.numeric(fit$mod$theta_hat)
			private$cached_values$hurdle_coefficients = fit$mod$hurdle_b
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		}
	)
)
