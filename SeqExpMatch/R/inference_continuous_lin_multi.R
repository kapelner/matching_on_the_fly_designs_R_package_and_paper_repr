#' Lin (2013) Covariate-Adjusted OLS Inference for Continuous Responses
#'
#' Fits the Lin (2013) covariate-adjusted linear estimator for continuous responses.
#' The working model includes an intercept, treatment
#' indicator, centered covariates, and treatment-by-centered-covariate
#' interactions. Inference uses HC2 heteroskedasticity-robust standard errors.
#'
#' @export
InferenceContinMultLin = R6::R6Class("InferenceContinMultLin",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize a Lin (2013) inference object for a completed sequential
		#' experiment with a continuous response.
		#'
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with a
		#'   continuous response.
		#'   randomization inference.
		#' @param verbose Whether to print progress messages.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 20, response_type = "continuous")
		#' for (t in 1:20) {
		#' 	x_t = data.frame(x1 = rnorm(1), x2 = rnorm(1))
		#' 	w_t = seq_des$add_subject_to_experiment_and_assign(x_t)
		#' 	seq_des$add_subject_response(t, x_t$x1 + 0.5 * w_t + rnorm(1))
		#' }
		#' seq_des_inf = InferenceContinMultLin$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "continuous")
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes Lin's covariate-adjusted estimate of the treatment effect.
		#'
		#' @return The estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval using the HC2 robust
		#' standard error.
		#'
		#' @param alpha The confidence level in the computed confidence interval is
		#'   1 - \code{alpha}. The default is 0.05.
		#'
		#' @return A confidence interval for the treatment effect.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a two-sided p-value for the treatment effect.
		#'
		#' @param delta The null treatment effect. Defaults to 0.
		#'
		#' @return The approximate p-value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		build_lin_design_matrix = function(){
			X = as.matrix(private$get_X())
			p = ncol(X)
			if (p == 0L){
				X_lin = cbind(1, private$w)
				colnames(X_lin) = c("(Intercept)", "treatment")
				return(X_lin)
			}

			if (is.null(colnames(X))){
				colnames(X) = paste0("x", seq_len(p))
			}

			Xc = scale(X, center = TRUE, scale = FALSE)
			Xc = as.matrix(Xc)
			colnames(Xc) = colnames(X)

			X_int = Xc * private$w
			colnames(X_int) = paste0("treatment:", colnames(X))

			X_lin = cbind(1, private$w, Xc, X_int)
			colnames(X_lin)[1:2] = c("(Intercept)", "treatment")
			X_lin
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(invisible(NULL))
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			X_full = private$build_lin_design_matrix()
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat

			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$cached_values$beta_hat_T = NA_real_
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}

			mod = stats::lm.fit(X_fit, private$y)
			coef_hat = as.numeric(mod$coefficients)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat))){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}

			post_fit = tryCatch(
				ols_hc2_post_fit_cpp(X_fit, as.numeric(private$y), coef_hat, j_treat),
				error = function(e) NULL
			)
			if (is.null(post_fit)){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}

			coef_names = colnames(X_fit)
			beta_hat = post_fit$beta_hat
			vcov_hc2 = post_fit$vcov
			ssq_hat = post_fit$ssq_hat
			if (!is.finite(beta_hat) || !is.finite(ssq_hat) || ssq_hat < 0){
				beta_hat = NA_real_
				ssq_hat = NA_real_
			}

			names(coef_hat) = coef_names
			colnames(vcov_hc2) = rownames(vcov_hc2) = coef_names
			std_err = post_fit$std_err
			names(std_err) = coef_names
			z_vals = post_fit$z_vals
			names(z_vals) = coef_names

			private$cached_values$beta_hat_T = beta_hat
			private$cached_values$s_beta_hat_T = if (is.finite(ssq_hat)) sqrt(ssq_hat) else NA_real_
			private$cached_values$is_z = TRUE
			private$cached_values$df = nrow(X_fit) - ncol(X_fit)
			private$cached_values$full_coefficients = coef_hat
			private$cached_values$full_vcov = vcov_hc2

			summary_table = cbind(
				Value = coef_hat,
				`Std. Error` = std_err,
				`z value` = z_vals,
				`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
			)
			private$cached_values$summary_table = summary_table
		}
	)
)
