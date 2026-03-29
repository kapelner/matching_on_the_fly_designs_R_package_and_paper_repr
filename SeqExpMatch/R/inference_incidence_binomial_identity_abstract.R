#' Binomial Identity-Link Inference for Binary Responses
#'
#' @description
#' Internal base class for binomial identity-link regression inference on
#' incidence outcomes. The treatment effect is reported on the risk-difference
#' scale.
#'
#' @keywords internal
#' @noRd
InferenceIncidBinomialIdentityAbstract = R6::R6Class("InferenceIncidBinomialIdentityAbstract",
	inherit = InferenceAsymp,
	public = list(
		# @description
		# Returns the treatment effect estimate (risk difference).
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		# @description
		# Computes an asymptotic confidence interval for the risk difference.
		# @param alpha The confidence level is 1 - alpha. Default 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Computes an asymptotic two-sided p-value for the treatment effect.
		# @param delta Null treatment effect. Default 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	private = list(
		build_design_matrix = function() stop(class(self)[1], " must implement build_design_matrix()."),

		fit_constrained_binomial = function(X_fit, j_treat){
			fast_identity_binomial_regression_with_var_cpp(X_fit, as.numeric(private$y), j = j_treat)
		},

		method_label = function() "Binomial identity-link regression",

		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)) x_names = paste0("x", seq_len(p))
			x_names
		},

		set_failed_fit_cache = function(){
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$is_z = TRUE
			private$cached_values$df = NA_real_
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				stop("Binomial identity-link regression: could not compute a finite standard error.")
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			X_full = private$build_design_matrix()
			if (is.null(dim(X_full))) X_full = matrix(X_full, ncol = 2L)
			colnames(X_full) = c("(Intercept)", "treatment", if (ncol(X_full) > 2L) private$get_covariate_names() else NULL)

			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			mod = tryCatch(
				private$fit_constrained_binomial(X_fit, j_treat),
				error = function(e) NULL
			)
			if (is.null(mod) || !isTRUE(mod$converged) || !is.finite(mod$ssq_b_j)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = as.numeric(mod$b[j_treat])
			if (estimate_only) return(invisible(NULL))
			private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			private$cached_values$is_z = TRUE
			private$cached_values$df = nrow(X_fit) - ncol(X_fit)
		}
	)
)
