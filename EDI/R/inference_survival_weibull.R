#' Weibull AFT Inference for Survival Responses
#'
#' Fits a Weibull Accelerated Failure Time (AFT) model for survival responses
#' using the treatment indicator and, optionally, all recorded covariates as
#' predictors. The treatment effect is reported on the log-time-ratio scale.
#'
#' @export
InferenceSurvivalWeibullRegr = R6::R6Class("InferenceSurvivalWeibullRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(
		#' @description
		#' Initialize a Weibull-regression inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_default Whether to use smart optimizer start values by default.
		#' @param optimization_alg Character scalar specifying the optimization algorithm. 
		#'   Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_default = FALSE, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
				assertFormula(model_formula, null.ok = TRUE)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_default = smart_default)
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

			res = fast_weibull_regression_cpp(
				y = private$y, dead = private$dead, X = X,
				start_params = private$get_fit_warm_start_for_length("params", ncol(X) + 1L),
				smart_start = private$smart_default,
				estimate_only = TRUE, optimization_alg = private$optimization_alg
			)
			if (is.null(res) || !isTRUE(res$converged) || !is.finite(res$coefficients[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(c(as.numeric(res$coefficients), as.numeric(res$log_sigma)), "params")
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
			dead = as.numeric(private$dead)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					res = tryCatch(
						fast_weibull_regression_cpp(
							y = y, dead = dead, X = X_fit,
							start_params = start %||% private$get_fit_warm_start_for_length("params", ncol(X_fit) + 1L),
							fixed_idx = j_treat, fixed_values = delta,
							smart_start = private$smart_default,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$coefficients), log_sigma = res$log_sigma, neg_loglik = as.numeric(res$neg_ll))
				},
				extract_start = function(fit){
					c(as.numeric(fit$b), as.numeric(fit$log_sigma))
				},
				score = function(fit){
					params = c(as.numeric(fit$b), as.numeric(fit$log_sigma))
					get_weibull_regression_score_cpp(X_fit, y, dead, params)
				},
				observed_information = function(fit){
					params = c(as.numeric(fit$b), as.numeric(fit$log_sigma))
					-get_weibull_regression_hessian_cpp(X_fit, y, dead, params)
				},
				fisher_information = function(fit){
					params = c(as.numeric(fit$b), as.numeric(fit$log_sigma))
					-get_weibull_regression_hessian_cpp(X_fit, y, dead, params)
				},
				information = function(fit){
					params = c(as.numeric(fit$b), as.numeric(fit$log_sigma))
					-get_weibull_regression_hessian_cpp(X_fit, y, dead, params)
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
					start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + 1L)
					res = fast_weibull_regression_cpp(
						y = private$y, dead = private$dead, X = X_fit,
						start_params = start_params,
						smart_start = private$smart_default,
						estimate_only = estimate_only, optimization_alg = private$optimization_alg
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(
						b = as.numeric(res$coefficients),
						log_sigma = as.numeric(res$log_sigma),
						neg_loglik = as.numeric(res$neg_ll),
						ssq_b_2 = if (estimate_only || is.null(res$vcov)) NA_real_ else {
							if (nrow(res$vcov) >= 2L) res$vcov[2L, 2L] else NA_real_
						}
					)
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
				}
			)

			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(c(as.numeric(attempt$fit$b), as.numeric(attempt$fit$log_sigma)), "params")
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
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
