#' Probit Regression Inference for Incidence Responses
#'
#' Fits a probit regression for binary (incidence) responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidProbitRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidProbitRegr = R6::R6Class("InferenceIncidProbitRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize a probit-regression inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values by default.
		#' @param optimization_alg  Optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = TRUE, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
				assertFormula(model_formula, null.ok = TRUE)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE, default = "newton_raphson")
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = as.numeric(private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights))
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit, keep){
					fit = tryCatch(
						stats::glm.fit(
							x = cbind(`(Intercept)` = 1, X_fit),
							y = as.numeric(private$y),
							weights = row_weights,
							family = stats::binomial(link = "probit")
						),
						error = function(e) NULL
					)
					if (is.null(fit) || is.null(fit$coefficients) || length(fit$coefficients) < 2L) return(NULL)
					list(
						b = as.numeric(fit$coefficients[-1L]),
						full_b = as.numeric(fit$coefficients),
						j_treat = match(1L, keep),
						params = as.numeric(fit$coefficients),
						fisher_information = tryCatch(solve(stats::vcov(stats::glm(y ~ .,
							data = data.frame(y = as.numeric(private$y), cbind(`(Intercept)` = 1, X_fit)[, -1, drop = FALSE]),
							weights = row_weights,
							family = stats::binomial(link = "probit")))), error = function(e) NULL)
					)
				},
				fit_ok = function(mod, X_fit, keep){
					!is.null(mod) && !is.na(mod$j_treat) && length(mod$b) >= mod$j_treat && is.finite(mod$b[mod$j_treat])
				}
			)
			if (is.null(attempt$fit) || is.null(attempt$fit$b) || length(attempt$fit$b) < 1L || !is.finite(attempt$fit$b[1L])) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			private$set_fit_warm_start(as.numeric(attempt$fit$params), "params", fisher = attempt$fit$fisher_information)
			private$cached_values$beta_hat_T = as.numeric(attempt$fit$b[attempt$fit$j_treat])
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		best_X_colnames = NULL,
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
				X = matrix(private$w, ncol = 1L)
				colnames(X) = "treatment"
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(treatment = private$w, X_cov)
			}
			n_params = ncol(X) + 1L
			ws_args = private$get_backend_warm_start_args(n_params)
			res = fast_ordinal_probit_regression_cpp(
				X = X,
				y = as.numeric(private$y),
				warm_start_params = ws_args$start_params,
				warm_start_fisher_info = ws_args$warm_start_fisher_info,
				smart_cold_start = private$smart_cold_start_default,
				optimization_alg = private$optimization_alg
			)
			if (is.null(res) || length(res$b) < 1L || !is.finite(res$b[1L])){
				return(NA_real_)
			}
			private$set_fit_warm_start(as.numeric(res$params), "params", fisher = res$fisher_information)
			as.numeric(res$b[1L])
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
			params_null = as.numeric(null_fit$params)
			p_reg  = ncol(spec$X)
			b_null = params_null[seq_len(p_reg)]
			kappa  = params_null[p_reg + 1L]
			eta    = as.numeric(spec$X %*% b_null)
			mu     = pmin(pmax(pnorm(eta - kappa), 0), 1)
			y_sim  = private$simulate_param_boot_bernoulli_y(mu)
			if (is.null(y_sim)) return(NULL)
			X_fit  = spec$X
			j      = spec$j
			
			# Parametric bootstrap: use observed fit as anchor
			ws_args = private$get_backend_warm_start_args(ncol(X_fit) + 1L)
			full_b = tryCatch(
				fast_ordinal_probit_regression_cpp(
					X = X_fit, y = y_sim,
					warm_start_params = ws_args$start_params,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = TRUE,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_b) || length(full_b) == 0L) return(NULL)
			full_fit_boot = list(params = as.numeric(full_b$params), neg_loglik = as.numeric(full_b$neg_loglik))
			if (!is.finite(full_fit_boot$neg_loglik)) return(NULL)
			list(
				worker_data = list(y = y_sim),
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					ws_args_null = private$get_backend_warm_start_args(ncol(X_fit) + 1L)
					res = tryCatch(
						fast_ordinal_probit_regression_cpp(
							X = X_fit, y = y_sim,
							warm_start_params = start %||% full_fit_boot$params,
							warm_start_fisher_info = ws_args_null$warm_start_fisher_info,
							smart_cold_start = TRUE,
							optimization_alg = private$optimization_alg,
							fixed_idx = j, fixed_values = d
						),
						error = function(e) NULL
					)
					if (is.null(res) || length(res) == 0L) return(NULL)
					list(params = as.numeric(res$params), neg_loglik = as.numeric(res$neg_loglik))
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik)
			)
		},
		supports_fisher_information = function(){
			TRUE
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			full_fit = list(params = ctx$full_params, neg_loglik = ctx$full_neg_loglik)
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = full_fit,
				fit_null = function(delta, start = NULL){
					ws_args = private$get_backend_warm_start_args(length(ctx$full_params))
					res = tryCatch(
						fast_ordinal_probit_regression_cpp(
							X_fit,
							y,
							warm_start_params = start %||% ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						),
						error = function(e) NULL
					)
					if (is.null(res) || length(res) == 0L) return(NULL)
					list(params = as.numeric(res$params), neg_loglik = as.numeric(res$neg_loglik), fisher_information = res$fisher_information)
				},
				extract_start = function(fit){
					as.numeric(fit$params)
				},
				score = function(fit){
					get_ordinal_probit_regression_score_cpp(X_fit, y, as.numeric(fit$params))
				},
				observed_information = function(fit){
					-get_ordinal_probit_regression_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				fisher_information = function(fit){
					-get_ordinal_probit_regression_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				information = function(fit){
					-get_ordinal_probit_regression_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik)
				}
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit, keep){
					j_treat = match(1L, keep)
					n_params = ncol(X_fit) + 1L
					ws_args = private$get_backend_warm_start_args(n_params)
					if (estimate_only) {
						res = fast_ordinal_probit_regression_cpp(
							X = X_fit,
							y = as.numeric(private$y),
							warm_start_params = ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default,
							optimization_alg = private$optimization_alg
						)
						if (is.null(res) || length(res) == 0L) return(NULL)
						list(
							b = as.numeric(res$b),
							ssq_b_j = NA_real_,
							j_treat = j_treat,
							params = as.numeric(res$params),
							neg_loglik = as.numeric(res$neg_loglik),
							fisher_information = res$fisher_information
						)
					} else {
						res = fast_ordinal_probit_regression_with_var_cpp(
							X = X_fit,
							y = as.numeric(private$y),
							warm_start_params = ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							smart_cold_start = private$smart_cold_start_default,
							optimization_alg = private$optimization_alg
						)
						if (is.null(res) || length(res$b) == 0L || is.na(res$b[1L])) return(NULL)
						list(
							b = as.numeric(res$b),
							ssq_b_j = as.numeric(res$ssq_b_j),
							j_treat = j_treat,
							params = as.numeric(res$params),
							neg_loglik = as.numeric(res$neg_loglik),
							fisher_information = res$fisher_information
						)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = mod$j_treat
					if (is.null(mod) || is.na(j_treat) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(as.numeric(attempt$fit$params), "params")
				private$best_X_colnames = setdiff(colnames(attempt$X), "treatment")
				n_alpha = length(attempt$fit$params) - ncol(attempt$X)
				j_treat = n_alpha + match(1L, attempt$keep)
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = as.integer(j_treat),
					full_params = as.numeric(attempt$fit$params),
					full_neg_loglik = as.numeric(attempt$fit$neg_loglik)
				)
				list(
					b = c(0, as.numeric(attempt$fit$b[attempt$fit$j_treat])),
					ssq_b_2 = as.numeric(attempt$fit$ssq_b_j)
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
				NULL
			}
		},
		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = matrix(private$w, ncol = 1L)
				colnames(X) = "treatment"
			} else {
				X = cbind(treatment = private$w, X_cov)
			}
			X
		}
	)
)
