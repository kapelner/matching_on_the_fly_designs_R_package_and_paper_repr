#' Proportional Odds Regression Inference for Ordinal Responses
#'
#' Fits a proportional odds regression for ordinal responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalPropOddsRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalPropOddsRegr = R6::R6Class("InferenceOrdinalPropOddsRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize a proportional-odds inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values by default.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
				assertFormula(model_formula, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
		,
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			X_fit = private$build_design_matrix()
			n_params = ncol(X_fit) + length(sort(unique(private$y))) - 1L
			ws_args = private$get_backend_warm_start_args(n_params)
			res = tryCatch(
				fast_ordinal_regression_weighted_cpp(
					X = X_fit,
					y = as.numeric(private$y),
					weights = as.numeric(row_weights),
					warm_start_params = ws_args$start_params,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)
			if (is.null(res) || length(res$b) < 1L || !is.finite(res$b[length(res$b)])){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(NA_real_)
			}
			private$set_fit_warm_start(res$params, "params", fisher = res$fisher_information)
			private$cached_values$beta_hat_T = as.numeric(res$b[length(res$b)])
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
				X = as.matrix(private$w)
				colnames(X) = "treatment"
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(treatment = private$w, X_cov)
			}
			n_params = ncol(X) + length(sort(unique(private$y))) - 1L
			ws_args = private$get_backend_warm_start_args(n_params)
			res = fast_ordinal_regression_cpp(
				X = X, y = as.numeric(private$y),
				warm_start_params = ws_args$start_params,
				warm_start_fisher_info = ws_args$warm_start_fisher_info,
				smart_cold_start = private$smart_cold_start_default
			)
			if (is.null(res) || length(res$b) < 1L || !is.finite(res$b[length(res$b)])){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$params, "params", fisher = res$fisher_information)
			as.numeric(res$b[length(res$b)])
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
			y_sim       = private$simulate_param_boot_ordinal_y(spec$X, params_null, spec$y, stats::plogis)
			if (is.null(y_sim)) return(NULL)
			X_fit    = spec$X
			j        = spec$j
			
			# Parametric bootstrap: use observed fit as anchor
			ws_args = private$get_backend_warm_start_args(length(params_null))
			full_res = tryCatch(
				fast_ordinal_regression_cpp(
					X_fit, y_sim,
					warm_start_params = ws_args$start_params,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = private$smart_cold_start_default
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || length(full_res$params) == 0L) return(NULL)
			full_fit_boot = list(params = as.numeric(full_res$params), neg_loglik = as.numeric(full_res$neg_loglik))
			if (!is.finite(full_fit_boot$neg_loglik)) return(NULL)
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					ws_args_null = private$get_backend_warm_start_args(length(params_null))
					res = tryCatch(
						fast_ordinal_regression_cpp(
							X_fit, y_sim,
							warm_start_params = start %||% full_fit_boot$params,
							warm_start_fisher_info = ws_args_null$warm_start_fisher_info,
							fixed_idx = j, fixed_values = d,
							smart_cold_start = TRUE
						),
						error = function(e) NULL
					)
					if (is.null(res) || length(res) == 0L) return(NULL)
					list(params = as.numeric(res$params), neg_loglik = as.numeric(res$neg_loglik))
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik)
			)
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
				X = X_fit, y = y, j = j_treat,
				full_fit = full_fit,
				fit_null = function(delta, start = NULL){
					ws_args = private$get_backend_warm_start_args(length(ctx$full_params))
					res = tryCatch(
						fast_ordinal_regression_cpp(
							X_fit, y,
							warm_start_params = start %||% ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							fixed_idx = j_treat, fixed_values = delta,
							smart_cold_start = private$smart_cold_start_default
						),
						error = function(e) NULL
					)
					if (is.null(res) || length(res) == 0) return(NULL)
					list(params = as.numeric(res$params), neg_loglik = as.numeric(res$neg_loglik), fisher_information = res$fisher_information)
				},
				extract_start = function(fit){
					as.numeric(fit$params)
				},
				score = function(fit){
					get_ordinal_regression_score_cpp(X_fit, y, as.numeric(fit$params))
				},
				observed_information = function(fit){
					-get_ordinal_regression_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				fisher_information = function(fit){
					-get_ordinal_regression_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				information = function(fit){
					-get_ordinal_regression_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					n_params = ncol(X_fit) + length(sort(unique(private$y))) - 1L
					ws_args = private$get_backend_warm_start_args(n_params)
					if (estimate_only) {
						res = tryCatch(
							fast_ordinal_regression_cpp(
								X_fit, private$y,
								warm_start_params = ws_args$start_params,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default
							),
							error = function(e) NULL
						)
						list(b = res$b, ssq_b_j = NA_real_, params = res$params, neg_loglik = res$neg_loglik, fisher_information = res$fisher_information)
					} else {
						res = tryCatch(
							fast_ordinal_regression_with_var_cpp(
								X_fit, private$y,
								warm_start_params = ws_args$start_params,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default
							),
							error = function(e) NULL
						)
						list(b = res$b, ssq_b_j = res$ssq_b_j, params = res$params, neg_loglik = res$neg_loglik, fisher_information = res$fisher_information)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = length(mod$b)
					if (is.null(mod) || j_treat < 1L || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(attempt$fit$params, "params", fisher = attempt$fit$fisher_information)
				private$best_X_colnames = setdiff(colnames(attempt$X), "treatment")
				n_alpha = length(attempt$fit$params) - ncol(attempt$X)
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = as.integer(n_alpha + 1L),
					full_params = as.numeric(attempt$fit$params),
					full_neg_loglik = as.numeric(attempt$fit$neg_loglik)
				)
				list(b = c(0, attempt$fit$b[length(attempt$fit$b)]), ssq_b_2 = attempt$fit$ssq_b_j)
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
