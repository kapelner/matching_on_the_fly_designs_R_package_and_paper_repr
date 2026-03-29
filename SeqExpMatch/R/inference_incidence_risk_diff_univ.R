#' Univariate Risk-Difference Regression for Binary Responses
#'
#' @description
#' Fits a direct risk-difference estimator for binary (incidence) responses
#' using a linear probability model with HC2
#' heteroskedasticity-robust variance. The treatment effect is reported on the
#' risk-difference scale.
#'
#' @details
#' This class provides a direct additive treatment-effect model for binary
#' outcomes. With only the treatment indicator in the design matrix, the point
#' estimate coincides with the empirical risk difference.
#'
#' @export
InferenceIncidUnivRiskDiff = R6::R6Class("InferenceIncidUnivRiskDiff",
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize a risk-difference inference object for a completed
		#' design with a binary response.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with an
		#'   incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and
		#'   randomization inference.
		#' @param verbose Whether to print progress messages.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 20, response_type = "incidence")
		#' for (i in 1:20) {
		#' 	x_i = data.frame(x1 = rnorm(1), x2 = rnorm(1))
		#' 	w_i = seq_des$add_subject_to_experiment_and_assign(x_i)
		#' 	p_i = plogis(-0.8 + 0.5 * w_i)
		#' 	seq_des$add_subject_response(i, rbinom(1, 1, p_i))
		#' }
		#' seq_des_inf = InferenceIncidUnivRiskDiff$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "incidence")
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the direct risk-difference estimate of the treatment effect.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval using the HC2 robust
		#' standard error.
		#' @param alpha The confidence level in the computed confidence interval is
		#'   1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a two-sided p-value for the treatment effect.
		#' @param delta The null treatment effect on the risk-difference scale.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		},

		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				stop("Risk-difference regression: could not compute a finite HC2 standard error.")
			}
		},

		set_failed_fit_cache = function(){
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$is_z = TRUE
			private$cached_values$df = NA_real_
			private$cached_values$full_coefficients = NULL
			private$cached_values$full_vcov = NULL
			private$cached_values$summary_table = NULL
		},

		fit_risk_diff_model = function(X_fit, j_treat){
			mod = stats::lm.fit(X_fit, private$y)
			coef_hat = as.numeric(mod$coefficients)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat))){
				return(NULL)
			}

			post_fit = tryCatch(
				ols_hc2_post_fit_cpp(X_fit, as.numeric(private$y), coef_hat, j_treat),
				error = function(e) NULL
			)
			if (is.null(post_fit)){
				return(NULL)
			}

			coef_names = colnames(X_fit)
			names(coef_hat) = coef_names
			vcov_hc2 = post_fit$vcov
			colnames(vcov_hc2) = rownames(vcov_hc2) = coef_names
			std_err = post_fit$std_err
			names(std_err) = coef_names
			z_vals = post_fit$z_vals
			names(z_vals) = coef_names

			list(
				beta_hat = post_fit$beta_hat,
				se = post_fit$se,
				coefficients = coef_hat,
				vcov = vcov_hc2,
				summary_table = cbind(
					Value = coef_hat,
					`Std. Error` = std_err,
					`z value` = z_vals,
					`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
				)
			)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			X_full = private$build_design_matrix()
			if (is.null(dim(X_full))){
				X_full = matrix(X_full, ncol = 2L)
			}
			colnames(X_full) = c("(Intercept)", "treatment", if (ncol(X_full) > 2L) private$get_covariate_names() else NULL)

			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			fit = private$fit_risk_diff_model(X_fit, j_treat)
			if (is.null(fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = fit$beta_hat
			if (estimate_only) return(invisible(NULL))
			private$cached_values$s_beta_hat_T = fit$se
			private$cached_values$is_z = TRUE
			private$cached_values$df = nrow(X_fit) - ncol(X_fit)
			private$cached_values$full_coefficients = fit$coefficients
			private$cached_values$full_vcov = fit$vcov
			private$cached_values$summary_table = fit$summary_table
		}
	)
)
