#' OLS Inference for Continuous Responses
#'
#' Fits an ordinary least squares regression for continuous responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceContinOLS$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceContinOLS = R6::R6Class("InferenceContinOLS",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(
		#' @description Initialize an OLS inference object.
		#' @param des_obj A completed \code{Design} object with a continuous response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param max_resample_attempts Maximum number of times a single bootstrap replicate
		#'   may be redrawn when the drawn sample fails validity screening. Default \code{50L}.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, max_resample_attempts = 50L){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "continuous")
				assertCount(max_resample_attempts, positive = TRUE)
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			
			private$max_resample_attempts = max_resample_attempts
		},
		#' @description Computes the OLS estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			X_full = private$build_design_matrix()
			keep = is.finite(row_weights) & row_weights > 0 & is.finite(private$y)
			if (!any(keep)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			fit = tryCatch(
				stats::lm.wfit(
					x = X_full[keep, , drop = FALSE],
					y = as.numeric(private$y[keep]),
					w = as.numeric(row_weights[keep])
				),
				error = function(e) NULL
			)
			coef_hat = if (!is.null(fit)) as.numeric(stats::coef(fit)) else numeric(0)
			if (length(coef_hat) < 2L || !is.finite(coef_hat[2L])) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			private$cached_values$beta_hat_T = coef_hat[2L]
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval for the treatment effect.
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value for the treatment effect.
		#' @param delta The null difference to test against. Default is zero.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	private = list(
		max_resample_attempts = NULL,
		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			X
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			X_full = private$build_design_matrix()
			
			fit = tryCatch(stats::lm.fit(X_full, private$y), error = function(e) NULL)
			if (is.null(fit) || !is.finite(stats::coef(fit)[2])){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
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
		}
	)
)
