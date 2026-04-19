#' OLS Inference for Continuous Responses
#'
#' Fits an ordinary least squares regression for continuous responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceContinOLS = R6::R6Class("InferenceContinOLS",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize an OLS inference object.
		#' @param des_obj A completed \code{Design} object with a continuous response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		#' @param max_resample_attempts Maximum number of times a single bootstrap replicate
		#'   may be redrawn when the drawn sample fails validity screening. Default \code{50L}.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE, max_resample_attempts = 50L){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "continuous")
				assertCount(max_resample_attempts, positive = TRUE)
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
			private$max_resample_attempts = max_resample_attempts
		},

		#' @description
		#' Computes the OLS estimate of the treatment effect.
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
		max_resample_attempts = NULL,
		include_covariates = NULL,

		build_design_matrix = function(){
			if (private$include_covariates) {
				private$create_design_matrix()
			} else {
				X = cbind(1, private$w)
				colnames(X) = c("(Intercept)", "treatment")
				X
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			X_full = private$build_design_matrix()
			
			fit = tryCatch(stats::lm.fit(X_full, private$y), error = function(e) NULL)
			if (is.null(fit) || !is.finite(stats::coef(fit)[2])){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = FALSE
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = as.numeric(stats::coef(fit)[2])
			if (estimate_only) return(invisible(NULL))

			# Standard OLS variance calculation
			res = stats::residuals(fit)
			rss = sum(res^2)
			df  = nrow(X_full) - ncol(X_full)
			sig2 = rss / df
			
			# (X'X)^-1
			v_cov = sig2 * chol2inv(fit$qr$qr[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE])
			
			private$cached_values$s_beta_hat_T = sqrt(v_cov[2, 2])
			private$cached_values$df = df
			private$cached_values$is_z = FALSE
		}
	)
)
