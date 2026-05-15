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
		#' @param smart_default Whether to use smart optimizer start values by default.
		#' @param optimization_alg  Optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_default = FALSE, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_default = smart_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	),
	private = list(
		best_X_colnames = NULL,
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
				fast_neg_bin_cpp(
					X = X, y = as.integer(private$y),
					warm_start_params = private$get_fit_warm_start_for_length("params", ncol(X) + 1L),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X) + 1L),
					smart_start = private$smart_default,
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
		supports_likelihood_tests = function(){
			TRUE
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
					res = tryCatch(
						fast_neg_bin_cpp(
							X = X_fit, y = y,
							warm_start_params = start %||% private$get_fit_warm_start_for_length("params", ncol(X_fit) + 1L),
							warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit) + 1L),
							fixed_idx = j_treat, fixed_values = delta,
							smart_start = private$smart_default,
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
					warm_start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + 1L)
					warm_fisher = private$get_fit_warm_start_fisher(ncol(X_fit) + 1L)
					if (estimate_only) {
						res = tryCatch(
							fast_neg_bin_cpp(
								X = X_fit, y = as.integer(private$y),
								warm_start_params = warm_start_params,
								warm_start_fisher_info = warm_fisher,
								smart_start = private$smart_default,
								optimization_alg = private$optimization_alg
							),
							error = function(e) NULL
						)
						if (is.null(res) || !isTRUE(res$converged)) return(NULL)
						list(b = as.numeric(res$b), ssq_b_j = NA_real_, j_treat = j_treat,
						     theta_hat = res$theta_hat, neg_loglik = -as.numeric(res$logLik),
						     fisher_information = res$fisher_information)
					} else {
						res = tryCatch(
							fast_neg_bin_with_var_cpp(
								X = X_fit, y = as.integer(private$y),
								warm_start_params = warm_start_params,
								warm_start_fisher_info = warm_fisher,
								smart_start = private$smart_default,
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
