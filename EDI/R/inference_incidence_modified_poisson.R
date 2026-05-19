#' Modified Poisson Regression Inference for Incidence Responses
#'
#' Fits a modified Poisson regression (Zou 2004) for binary (incidence) responses
#' using the treatment indicator and, optionally, all recorded covariates as
#' predictors. This model provides an alternative to log-binomial regression for
#' estimating risk ratios.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidModifiedPoisson$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidModifiedPoisson = R6::R6Class("InferenceIncidModifiedPoisson",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
				
		#' @description Initialize a modified Poisson regression inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Compute the treatment effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Row weights for the bootstrap sample.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = private$build_design_matrix(),
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					res = fast_poisson_regression_weighted_cpp(
						X = X_fit,
						y = private$y,
						weights = row_weights,
						warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
						smart_cold_start = private$smart_cold_start_default,
						warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit))
					)
					list(b = res$b, ssq_b_j = NA_real_, j_treat = j_treat, mod = res, XtWX = res$XtWX)
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = mod$j_treat
					!is.null(mod) && length(mod$b) >= j_treat && is.finite(mod$b[j_treat])
				}
			)
			private$cached_mod = attempt$fit$mod %||% attempt$fit
			if (is.null(attempt$fit) || is.null(attempt$fit$b)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				return(NA_real_)
			}
			j_treat = attempt$fit$j_treat %||% 2L
			private$cached_values$beta_hat_T = as.numeric(attempt$fit$b[j_treat])
			private$cached_values$s_beta_hat_T = NA_real_
			private$set_fit_warm_start(
				as.numeric(attempt$fit$b),
				"beta",
				fisher = attempt$fit$XtWX %||% attempt$fit$fisher_information
			)
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		best_X_colnames = NULL,
		cached_mod = NULL,
		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			X
		},
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_cols = private$best_X_colnames
			X_data = private$get_X()
			if (length(X_cols) == 0L){
				X = cbind(1, private$w)
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(1, treatment = private$w, X_cov)
			}
			res = tryCatch(
				fast_poisson_regression_cpp(
					X = X, y = as.numeric(private$y),
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X)),
					smart_cold_start = private$smart_cold_start_default,
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X))
				),
				error = function(e) NULL
			)

			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$b, "beta", fisher = res$XtWX)
			as.numeric(res$b[2])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		supports_lik_ratio_param_bootstrap = function(){
			TRUE
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null     = as.numeric(null_fit$b)
			mu         = pmax(exp(as.numeric(spec$X %*% b_null)), 0)
			y_sim      = as.numeric(rpois(length(mu), mu))
			X_fit      = spec$X
			j          = spec$j
			full_fit_b = tryCatch(
				fast_poisson_regression_cpp(
					X = X_fit, y = y_sim,
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)
			if (is.null(full_fit_b) || length(full_fit_b$b) < j || !is.finite(full_fit_b$b[j])) return(NULL)
			list(
				worker_data = list(y = y_sim),
				full_fit = full_fit_b,
				fit_null = function(d, start = NULL){
					tryCatch(
						fast_poisson_regression_with_var_cpp(
							X = X_fit, y = y_sim, j = j,
							warm_start_beta = start %||% full_fit_b$b,
							fixed_idx = j, fixed_values = d,
							smart_cold_start = TRUE
						),
						error = function(e) NULL
					)
				},
				neg_loglik = function(fit){
					eta_f = as.numeric(X_fit %*% as.numeric(fit$b))
					-sum(y_sim * eta_f - exp(eta_f) - lgamma(y_sim + 1))
				}
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					fast_poisson_regression_with_var_cpp(
						X = X_fit,
						y = y,
						j = j_treat,
						warm_start_beta = start %||% private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
						warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit)),
						fixed_idx = j_treat,
						fixed_values = delta,
						smart_cold_start = private$smart_cold_start_default
					)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					as.numeric(fit$score %||% get_poisson_regression_score_cpp(X_fit, y, as.numeric(fit$b)))
				},
				observed_information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				neg_loglik = function(fit){
					eta = as.numeric(X_fit %*% as.numeric(fit$b))
					-sum(y * eta - exp(eta) - lgamma(y + 1))
				}
			)
		},
		generate_mod = function(estimate_only = FALSE){
			# Use the common GLM fitting pattern
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = private$build_design_matrix(),
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit))
					warm_fisher = private$get_fit_warm_start_fisher(ncol(X_fit))
					if (estimate_only) {
						res = fast_poisson_regression_cpp(
							X = X_fit, y = private$y,
							warm_start_beta = warm_start_beta,
							smart_cold_start = private$smart_cold_start_default,
							warm_start_fisher_info = warm_fisher
						)
						list(b = res$b, ssq_b_j = NA_real_, j_treat = j_treat, mod = res, XtWX = res$XtWX)
					} else {
						res = fast_poisson_regression_with_var_cpp(
							X = X_fit, y = private$y, j = j_treat,
							warm_start_beta = warm_start_beta,
							smart_cold_start = private$smart_cold_start_default,
							warm_start_fisher_info = warm_fisher
						)
						res$j_treat = j_treat
						res
					}
				},

				fit_ok = function(mod, X_fit, keep){
					j_treat = mod$j_treat
					if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
				if (!is.null(attempt$fit)){
					private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
					private$cached_mod = attempt$fit$mod %||% attempt$fit
					private$cached_values$likelihood_test_context = list(X = attempt$X, j_treat = which(attempt$keep == 2L))
					if (!is.null(attempt$fit$b)) {
						private$set_fit_warm_start(attempt$fit$b, "beta", fisher = attempt$fit$XtWX %||% attempt$fit$fisher_information)
					}
				}
				attempt$fit
			}
	)
)

#' Multi-subject Modified Poisson Inference for Incidence Responses
#'
#' Historical public alias for the modified Poisson implementation.
#'
#' @export
