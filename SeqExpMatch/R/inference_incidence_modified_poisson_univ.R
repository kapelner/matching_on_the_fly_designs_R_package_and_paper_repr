#' Univariate Modified Poisson Inference for Binary Responses
#'
#' @description
#' Fits the classic modified Poisson model for binary (incidence) responses using
#' a Poisson log-link working model with Huber-White
#' sandwich variance. The treatment effect is reported on the log-risk-ratio
#' scale.
#'
#' @details
#' This is the practical risk-ratio approach popularized by Zou (2004). The
#' working model is Poisson with log link, but inference uses a robust sandwich
#' covariance rather than the naive Poisson model-based variance.
#'
#' @export
InferenceIncidUnivModifiedPoisson = R6::R6Class("InferenceIncidUnivModifiedPoisson",
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize a modified-Poisson inference object for a completed
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
		#' 	p_i = pmin(0.95, exp(-1 + 0.4 * w_i))
		#' 	seq_des$add_subject_response(i, rbinom(1, 1, p_i))
		#' }
		#' seq_des_inf = InferenceIncidUnivModifiedPoisson$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "incidence")
			super$initialize(des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the modified-Poisson estimate of the treatment effect on the
		#' log-risk-ratio scale.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval using the sandwich
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
		#' @param delta The null treatment effect on the log-risk-ratio scale.
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
				stop("Modified Poisson estimator: could not compute a finite sandwich standard error.")
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

		fit_modified_poisson = function(X_fit, j_treat){
			mod = tryCatch(
				suppressWarnings(stats::glm.fit(x = X_fit, y = as.numeric(private$y), family = stats::poisson(link = "log"))),
				error = function(e) NULL
			)
			if (is.null(mod) || !isTRUE(mod$converged)){
				return(NULL)
			}

			coef_hat = as.numeric(mod$coefficients)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat))){
				return(NULL)
			}

			mu_hat = as.numeric(mod$fitted.values)
			if (length(mu_hat) != nrow(X_fit) || any(!is.finite(mu_hat)) || any(mu_hat <= 0)){
				return(NULL)
			}

			post_fit = tryCatch(
				glm_sandwich_post_fit_cpp(
					X_fit = X_fit,
					y = as.numeric(private$y),
					coef_hat = coef_hat,
					mu_hat = mu_hat,
					working_weights = mu_hat,
					j_treat = j_treat
				),
				error = function(e) NULL
			)
			if (is.null(post_fit)){
				return(NULL)
			}

			coef_names = colnames(X_fit)
			names(coef_hat) = coef_names
			vcov_robust = post_fit$vcov
			colnames(vcov_robust) = rownames(vcov_robust) = coef_names
			std_err = post_fit$std_err
			names(std_err) = coef_names
			z_vals = post_fit$z_vals
			names(z_vals) = coef_names

			list(
				beta_hat = post_fit$beta_hat,
				se = post_fit$se,
				coefficients = coef_hat,
				vcov = vcov_robust,
				summary_table = cbind(
					Value = coef_hat,
					`Std. Error` = std_err,
					`z value` = z_vals,
					`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
				)
			)
		},

		shared = function(){
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

			fit = private$fit_modified_poisson(X_fit, j_treat)
			if (is.null(fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = fit$beta_hat
			private$cached_values$s_beta_hat_T = fit$se
			private$cached_values$is_z = TRUE
			private$cached_values$df = nrow(X_fit) - ncol(X_fit)
			private$cached_values$full_coefficients = fit$coefficients
			private$cached_values$full_vcov = fit$vcov
			private$cached_values$summary_table = fit$summary_table
		}
	)
)
