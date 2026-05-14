#' Beta Regression Inference for Proportion Responses
#'
#' Fits a beta regression for proportion responses (constrained to (0, 1)) using
#' the treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'proportion')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferencePropBetaRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferencePropBetaRegr = R6::R6Class("InferencePropBetaRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize a beta-regression inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg Character scalar specifying the optimization algorithm. 
		#'   Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "proportion")
				assertFormula(model_formula, null.ok = TRUE)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
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
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			n_params = ncol(X) + 1L
			res = fast_beta_regression_cpp(
				X = X, y = as.numeric(private$y),
				start_beta = private$get_fit_warm_start_for_length("beta", ncol(X)),
				warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
				optimization_alg = private$optimization_alg
			)
			if (is.null(res) || !is.finite(res$coefficients[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(c(as.numeric(res$coefficients), log(as.numeric(res$phi))), "params", fisher = res$fisher_information)
			as.numeric(res$coefficients[2])
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
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					res = fast_beta_regression_cpp(
						X_fit, y,
						fixed_idx = j_treat, fixed_values = delta,
						start_beta = start[1:ncol(X_fit)],
						warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit) + 1L),
						optimization_alg = private$optimization_alg
					)
					if (is.null(res)) return(NULL)
					list(b = as.numeric(res$coefficients), phi = res$phi, neg_loglik = res$neg_loglik, fisher_information = res$fisher_information)
				},
				extract_start = function(fit){
					c(as.numeric(fit$b), log(as.numeric(fit$phi)))
				},
				score = function(fit){
					params = c(as.numeric(fit$b), log(as.numeric(fit$phi)))
					get_beta_regression_score_cpp(X_fit, y, params)
				},
				observed_information = function(fit){
					params = c(as.numeric(fit$b), log(as.numeric(fit$phi)))
					-get_beta_regression_hessian_cpp(X_fit, y, params)
				},
				fisher_information = function(fit){
					params = c(as.numeric(fit$b), log(as.numeric(fit$phi)))
					-get_beta_regression_hessian_cpp(X_fit, y, params)
				},
				information = function(fit){
					params = c(as.numeric(fit$b), log(as.numeric(fit$phi)))
					-get_beta_regression_hessian_cpp(X_fit, y, params)
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit){
					n_params = ncol(X_fit) + 1L
					start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit))
					warm_fisher = private$get_fit_warm_start_fisher(n_params)
					if (estimate_only) {
						res = fast_beta_regression_cpp(
							X_fit, private$y,
							start_beta = start_beta,
							warm_start_fisher_info = warm_fisher,
							optimization_alg = private$optimization_alg
						)
						if (is.null(res)) return(NULL)
						list(b = res$coefficients, ssq_b_2 = NA_real_, phi = res$phi, neg_loglik = res$neg_loglik, fisher_information = res$fisher_information)
					} else {
						res = fast_beta_regression_with_var_cpp(
							X_fit, private$y,
							start_beta = start_beta,
							warm_start_fisher_info = warm_fisher,
							optimization_alg = private$optimization_alg
						)
						if (is.null(res)) return(NULL)
						list(b = res$coefficients,
						     ssq_b_2 = if (!is.null(res$vcov) && nrow(res$vcov) >= 2L) res$vcov[2L, 2L] else NA_real_,
						     phi = res$phi, neg_loglik = res$neg_loglik, fisher_information = res$fisher_information)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$set_fit_warm_start(c(as.numeric(attempt$fit$b), log(as.numeric(attempt$fit$phi))), "params", fisher = attempt$fit$fisher_information)
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = which(attempt$keep == 2L)
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		},
		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			X
		}
	)
)
