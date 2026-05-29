#' Log-Binomial Regression Inference for Incidence Responses
#'
#' Fits a log-binomial regression for binary (incidence) responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidLogBinomial$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidLogBinomial = R6::R6Class("InferenceIncidLogBinomial",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(

		#' @description Initialize a log-binomial regression inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose               Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			X_data = private$get_X()
			X = if (is.null(X_data) || ncol(X_data) == 0) {
				cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				cbind(`(Intercept)` = 1, treatment = private$w, X_data)
			}
			res = tryCatch(
				fast_log_binomial_regression_weighted_cpp(
					X = X,
					y = as.numeric(private$y),
					weights = as.numeric(row_weights),
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X)),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X)),
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)
			if (is.null(res) || length(res$b) < 2L || !is.finite(res$b[2L])){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			private$set_fit_warm_start(res$b, "beta", fisher = res$fisher_information, force_pd = TRUE)
			private$cached_values$beta_hat_T = as.numeric(res$b[2L])
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		best_X_colnames = NULL,
		logbin_X_full_cache = NULL,
		logbin_w_cache = NULL,
		get_complexity_tier = function() "heavy",
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
			ws_args = private$get_backend_warm_start_args(ncol(X))
			res = tryCatch(
				fast_log_binomial_regression_cpp(
					X = X, y = as.numeric(private$y),
					warm_start_beta = ws_args$warm_start_beta,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)

			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$b, "beta", fisher = res$fisher_information, force_pd = TRUE)
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
			mu         = pmin(pmax(exp(as.numeric(spec$X %*% b_null)), 0), 1)
			y_sim      = as.numeric(rbinom(length(mu), 1L, mu))
			X_fit      = spec$X
			j          = spec$j

			# Parametric bootstrap: use observed fit as anchor
			ws_args = private$get_backend_warm_start_args(ncol(X_fit))
			full_fit_b = tryCatch(
				fast_log_binomial_regression_cpp(
					X = X_fit, y = y_sim,
					warm_start_beta = ws_args$warm_start_beta,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)
			if (is.null(full_fit_b) || !isTRUE(full_fit_b$converged)) return(NULL)
			if (length(full_fit_b$b) < j || !is.finite(full_fit_b$b[j])) return(NULL)
			list(
				full_fit = full_fit_b,
				fit_null = function(d, start = NULL){
					ws_args_null = private$get_backend_warm_start_args(ncol(X_fit))
					res = tryCatch(
						fast_log_binomial_regression_cpp(
							X = X_fit, y = y_sim,
							warm_start_beta = start %||% full_fit_b$b,
							warm_start_fisher_info = ws_args_null$warm_start_fisher_info,
							fixed_idx = j, fixed_values = d,
							smart_cold_start = TRUE
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					res
				},
				neg_loglik = function(fit){
					eta_f  = as.numeric(X_fit %*% as.numeric(fit$b))
					mu_fit = exp(eta_f)
					-sum(y_sim * log(pmax(mu_fit, 1e-15)) + (1 - y_sim) * log(pmax(1 - mu_fit, 1e-15)))
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
				X = X_fit, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					ws_args = private$get_backend_warm_start_args(ncol(X_fit))
					res = tryCatch(
						fast_log_binomial_regression_cpp(
							X_fit, y,
							warm_start_beta = start %||% ws_args$warm_start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							fixed_idx = j_treat, fixed_values = delta,
							smart_cold_start = private$smart_cold_start_default
						),
						error = function(e) NULL
					)

					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					res
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					get_log_binomial_regression_score_cpp(X_fit, y, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_log_binomial_regression_hessian_cpp(X_fit, y, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_log_binomial_regression_hessian_cpp(X_fit, y, as.numeric(fit$b))
				},
				information = function(fit){
					-get_log_binomial_regression_hessian_cpp(X_fit, y, as.numeric(fit$b))
				},
				neg_loglik = function(fit){
					eta = as.numeric(X_fit %*% as.numeric(fit$b))
					mu = exp(eta)
					-sum(y * log(pmax(mu, 1e-15)) + (1 - y) * log(pmax(1 - mu, 1e-15)))
				}
			)
		},
		generate_mod = function(estimate_only = FALSE){
			if (is.null(private$logbin_X_full_cache) || !identical(private$w, private$logbin_w_cache)) {
				X_data = private$get_X()
				private$logbin_X_full_cache = if (is.null(X_data) || ncol(X_data) == 0) {
					cbind(`(Intercept)` = 1, treatment = private$w)
				} else {
					cbind(`(Intercept)` = 1, treatment = private$w, X_data)
				}
				private$logbin_w_cache = private$w
			}
			X_full = private$logbin_X_full_cache
			
			if (!private$harden) {
				ws_args = private$get_backend_warm_start_args(ncol(X_full))
				if (estimate_only) {
					res = tryCatch(
						fast_log_binomial_regression_cpp(
							X_full, private$y,
							warm_start_beta = ws_args$warm_start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default,
							estimate_only = TRUE
						),
						error = function(e) NULL
					)
					if (is.null(res)) return(NULL)
					res$beta_hat_T = as.numeric(res$b[2L])
					res$ssq_b_j = NA_real_
					res$ssq_b_2 = NA_real_
				} else {
					res = tryCatch(
						fast_log_binomial_regression_with_var_cpp(
							X_full, private$y, j = 2L,
							warm_start_beta = ws_args$warm_start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default
						),
						error = function(e) NULL
					)
					if (is.null(res)) return(NULL)
					res$j_treat = 2L
					res$beta_hat_T = as.numeric(res$b[2L])
					res$ssq_b_2 = res$ssq_b_j
				}
				private$best_X_colnames = setdiff(colnames(X_full), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(
					X = X_full,
					j_treat = 2L,
					full_neg_loglik = res$neg_ll
				)
				return(res)
			}

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L, # intercept and treatment
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					ws_args = private$get_backend_warm_start_args(ncol(X_fit))
					if (estimate_only) {
						res = tryCatch(
							fast_log_binomial_regression_cpp(
								X = X_fit, y = private$y,
								warm_start_beta = ws_args$warm_start_beta,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default
							),
							error = function(e) NULL
						)
						if (is.null(res)) return(NULL)
						list(b = res$b, ssq_b_j = NA_real_, j_treat = j_treat, fisher_information = res$fisher_information, neg_ll = res$neg_ll)
					} else {
						res = tryCatch(
							fast_log_binomial_regression_with_var_cpp(
								X = X_fit, y = private$y, j = j_treat,
								warm_start_beta = ws_args$warm_start_beta,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default
							),
							error = function(e) NULL
						)
						if (is.null(res)) return(NULL)
						res$j_treat = j_treat
						res$ssq_b_2 = res$ssq_b_j
						res
					}
				},

				fit_ok = function(mod, X_fit, keep){
					j_treat = mod$j_treat
					if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j %||% mod$ssq_b_2)
				}
			)
			if (!is.null(attempt$fit)){
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = attempt$fit$j_treat,
					full_neg_loglik = attempt$fit$neg_ll
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		}
	)
)
