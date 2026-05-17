#' Negative Binomial Regression Inference for Count Responses
#'
#' Fits a negative binomial regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountNegBin$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountNegBin = R6::R6Class("InferenceCountNegBin",
	lock_objects = FALSE,
	inherit = InferenceCountLikelihood,
	public = list(
				
		#' @description Initialize a negative binomial regression inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart optimizer start values by default.
		#' @param optimization_alg  Optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = TRUE, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
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
				X = cbind(1, private$w)
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(1, treatment = private$w, X_cov)
			}
			ws_args = private$get_backend_warm_start_args(ncol(X) + 1L)
			res = tryCatch(
				fast_neg_bin_cpp(
					X = X, y = as.integer(private$y),
					warm_start_params = ws_args$start_params,
					warm_start_fisher_info = ws_args$warm_start_fisher_info,
					smart_cold_start = private$smart_cold_start_default,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(c(as.numeric(res$b), log(as.numeric(res$theta_hat))), "params", fisher = res$fisher_information)
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
			b_null   = as.numeric(null_fit$b)
			theta    = as.numeric(null_fit$theta_hat)
			if (!is.finite(theta) || theta <= 0) return(NULL)
			mu       = pmax(exp(as.numeric(spec$X %*% b_null)), 0)
			y_sim    = as.integer(rnbinom(length(mu), size = theta, mu = mu))
			X_fit    = spec$X
			j        = spec$j
			full_res = tryCatch(
				fast_neg_bin_cpp(
					X = X_fit, y = y_sim,
					smart_cold_start = private$smart_cold_start_default,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged) || length(full_res$b) < j || !is.finite(full_res$b[j])) return(NULL)
			full_fit_boot = list(
				b          = as.numeric(full_res$b),
				theta_hat  = full_res$theta_hat,
				neg_loglik = -as.numeric(full_res$logLik)
			)
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					res = tryCatch(
						fast_neg_bin_cpp(
							X = X_fit, y = y_sim,
							warm_start_params = start %||% c(as.numeric(full_res$b), log(as.numeric(full_res$theta_hat))),
							fixed_idx = j, fixed_values = d,
							smart_cold_start = TRUE,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$b), theta_hat = res$theta_hat, neg_loglik = -as.numeric(res$logLik))
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik)
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.integer(private$y)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					ws_args = private$get_backend_warm_start_args(ncol(X_fit) + 1L)
					res = tryCatch(
						fast_neg_bin_cpp(
							X = X_fit, y = y,
							warm_start_params = start %||% ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info,
							fixed_idx = j_treat, fixed_values = delta,
							smart_cold_start = private$smart_cold_start_default,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$b), theta_hat = res$theta_hat, neg_loglik = -as.numeric(res$logLik), fisher_information = res$fisher_information)
				},
				extract_start = function(fit){
					c(as.numeric(fit$b), log(as.numeric(fit$theta_hat)))
				},
				score = function(fit){
					params = c(as.numeric(fit$b), log(as.numeric(fit$theta_hat)))
					get_negbin_regression_score_cpp(X_fit, y, params)
				},
				observed_information = function(fit){
					params = c(as.numeric(fit$b), log(as.numeric(fit$theta_hat)))
					-get_negbin_regression_hessian_cpp(X_fit, y, params)
				},
				fisher_information = function(fit){
					params = c(as.numeric(fit$b), log(as.numeric(fit$theta_hat)))
					-get_negbin_regression_hessian_cpp(X_fit, y, params)
				},
				information = function(fit){
					params = c(as.numeric(fit$b), log(as.numeric(fit$theta_hat)))
					-get_negbin_regression_hessian_cpp(X_fit, y, params)
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_data = private$get_X()
			X_full = if (is.null(X_data) || ncol(X_data) == 0) {
				cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				cbind(`(Intercept)` = 1, treatment = private$w, X_data)
			}
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					ws_args = private$get_backend_warm_start_args(ncol(X_fit) + 1L)
					if (estimate_only) {
						res = tryCatch(
							fast_neg_bin_cpp(
								X = X_fit, y = as.integer(private$y),
								warm_start_params = ws_args$start_params,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default,
								optimization_alg = private$optimization_alg
							),
							error = function(e) NULL
						)
						if (is.null(res) || !isTRUE(res$converged)) return(NULL)
						list(b = as.numeric(res$b), ssq_b_2 = NA_real_, j_treat = j_treat,
						     theta_hat = res$theta_hat, neg_loglik = -as.numeric(res$logLik),
						     fisher_information = res$fisher_information)
					} else {
						res = tryCatch(
							fast_neg_bin_with_var_cpp(
								X = X_fit, y = as.integer(private$y),
								warm_start_params = ws_args$start_params,
								warm_start_fisher_info = ws_args$warm_start_fisher_info,
								smart_cold_start = private$smart_cold_start_default,
								optimization_alg = private$optimization_alg
							),
							error = function(e) NULL
						)
						if (is.null(res) || !isTRUE(res$converged)) return(NULL)
						hess = res$hess_fisher_info_matrix
						vcov = tryCatch(solve(hess), error = function(e) matrix(NA_real_, nrow(hess), ncol(hess)))
						ssq_b_j = if (j_treat <= ncol(X_fit)) as.numeric(vcov[j_treat, j_treat]) else NA_real_
						list(b = as.numeric(res$b), ssq_b_j = ssq_b_j, j_treat = j_treat,
						     theta_hat = res$theta_hat, neg_loglik = -as.numeric(res$logLik),
						     fisher_information = res$hess_fisher_info_matrix)
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
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = attempt$fit$j_treat
				)
				# Important: pack parameters for the base class shared logic
				attempt$fit$params = c(as.numeric(attempt$fit$b), log(as.numeric(attempt$fit$theta_hat)))
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		}
	)
)
