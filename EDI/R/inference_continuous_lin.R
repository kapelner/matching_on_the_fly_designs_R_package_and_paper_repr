#' Lin (2013) Covariate-Adjusted OLS Inference for Continuous Responses
#'
#' Fits the Lin (2013) covariate-adjusted linear estimator for continuous responses.
#' The working model includes an intercept, treatment indicator and, optionally,
#' centered covariates and treatment-by-centered-covariate interactions.
#' Inference uses HC2 heteroskedasticity-robust standard errors.
#'
#' @export
InferenceContinLin = R6::R6Class("InferenceContinLin",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize a Lin (2013) inference object.
		#' @param des_obj A completed \code{Design} object with a continuous response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors with interactions. If \code{FALSE}, only the
		#'   treatment indicator is used (this reduces to a simple OLS). If \code{NULL}
		#'   (default), it is set to \code{TRUE} if the design contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "continuous")
				assertFlag(include_covariates, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			
			if (is.null(include_covariates)) {
				include_covariates = des_obj$has_covariates()
			}
			private$include_covariates = include_covariates
		},

		#' @description
		#' Computes Lin's estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval using HC2 robust standard error.
		#' @param alpha The confidence level. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a two-sided p-value for the treatment effect.
		#' @param delta The null treatment effect. Defaults to 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		include_covariates = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},

		build_lin_design_matrix = function(){
			if (!private$include_covariates) {
				X_lin = cbind(1, private$w)
				colnames(X_lin) = c("(Intercept)", "treatment")
				return(X_lin)
			}

			Xc_info = private$get_centered_covariates()
			if (is.null(Xc_info)){
				X_lin = cbind(1, private$w)
				colnames(X_lin) = c("(Intercept)", "treatment")
				return(X_lin)
			}

			Xc = Xc_info$Xc
			X_int = Xc * private$w
			colnames(X_int) = paste0("treatment:", colnames(Xc))

			X_lin = cbind(1, private$w, Xc, X_int)
			colnames(X_lin)[1:2] = c("(Intercept)", "treatment")
			X_lin
		},

		get_centered_covariates = function(){
			des_priv = private$des_obj_priv_int
			# Design-level cache: NULL = not computed, list() = computed with p=0, list(Xc=...) = computed
			if (!is.null(des_priv$lin_centered_covariates)) {
				cached = des_priv$lin_centered_covariates
				return(if (length(cached) == 0L) NULL else cached)
			}

			X = as.matrix(private$get_X())
			p = ncol(X)
			if (p == 0L){
				des_priv$lin_centered_covariates = list()  # sentinel: computed, no covariates
				return(NULL)
			}

			if (is.null(colnames(X))){
				colnames(X) = paste0("x", seq_len(p))
			}

			Xc = scale(X, center = TRUE, scale = FALSE)
			Xc = as.matrix(Xc)
			colnames(Xc) = colnames(X)
			result = list(Xc = Xc)
			des_priv$lin_centered_covariates = result
			result
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(invisible(NULL))
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && isTRUE(private$cached_values$lin_estimate_only_complete)) return(invisible(NULL))
			if (!estimate_only && isTRUE(private$cached_values$lin_full_complete)) return(invisible(NULL))

			X_full = private$build_lin_design_matrix()
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat

			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$cache_nonestimable_estimate("linear_model_design_unusable")
				if (estimate_only) private$cached_values$lin_estimate_only_complete = TRUE
				private$cached_values$lin_full_complete = TRUE
				return(invisible(NULL))
			}

			mod = stats::lm.fit(X_fit, private$y)
			coef_hat = as.numeric(mod$coefficients)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat))){
				private$cache_nonestimable_estimate("linear_model_coefficients_unavailable")
				if (estimate_only){
					private$cached_values$lin_estimate_only_complete = TRUE
				} else {
					private$cached_values$lin_full_complete = TRUE
				}
				return(invisible(NULL))
			}

			beta_hat = coef_hat[j_treat]
			private$cached_values$beta_hat_T = beta_hat
			if (estimate_only){
				private$cached_values$lin_estimate_only_complete = TRUE
				return(invisible(NULL))
			}

			post_fit = tryCatch(
				ols_hc2_post_fit_cpp(X_fit, as.numeric(private$y), coef_hat, j_treat),
				error = function(e) NULL
			)
			if (is.null(post_fit)){
				private$cache_nonestimable_estimate("linear_model_post_fit_unavailable")
				private$cached_values$lin_full_complete = TRUE
				return(invisible(NULL))
			}

			coef_names = colnames(X_fit)
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
			private$cached_values$lin_estimate_only_complete = TRUE

			summary_table = cbind(
				Value = coef_hat,
				`Std. Error` = std_err,
				`z value` = z_vals,
				`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
			)
			private$cached_values$summary_table = summary_table
			private$cached_values$lin_full_complete = TRUE
		}
	)
)
