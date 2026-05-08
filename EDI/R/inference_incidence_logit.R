#' Logistic Regression Inference for Incidence Responses
#'
#' Fits a logistic regression for binary (incidence) responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidLogRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidLogRegr = R6::R6Class("InferenceIncidLogRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(
				
		#' @description
		#' Initialize a logistic-regression inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_default Whether to use smart optimizer start values by default.
		#' @param optimization_alg  Optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_default = TRUE, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
				assertFormula(model_formula, null.ok = TRUE)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = TRUE, default = "lbfgs")
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_default = smart_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	),

	private = list(
		best_X_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}

			# Use the same design matrix structure as the original fit
			X_cols = private$best_X_colnames
			X_data = private$get_X()

			if (length(X_cols) == 0L){
				# Univariate case
				X = cbind(1, private$w)
			} else {
				# Multivariate case
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(1, treatment = private$w, X_cov)
			}

			res = fast_logistic_regression_cpp(
				X = X, y = as.numeric(private$y),
				start_beta = private$get_fit_warm_start_for_length("beta", ncol(X)),
				smart_start = private$smart_default,
				optimization_alg = private$optimization_alg
			)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$b, "beta")
			as.numeric(res$b[2])
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		supports_likelihood_tests = function(){
			TRUE
		},

		supports_fisher_information = function(){
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
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					fast_logistic_regression_with_var_cpp(
						X = X_fit,
						y = y,
						j = j_treat,
						start_beta = start %||% private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
						fixed_idx = j_treat,
						fixed_values = delta,
						smart_start = private$smart_default,
						optimization_alg = private$optimization_alg
					)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					get_logistic_regression_score_cpp(X_fit, y, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_logistic_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_logistic_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				information = function(fit){
					-get_logistic_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				neg_loglik = function(fit){
					eta = as.numeric(X_fit %*% as.numeric(fit$b))
					log_denom = ifelse(eta > 0, eta + log1p(exp(-eta)), log1p(exp(eta)))
					-sum(y * eta - log_denom)
				}
			)
		},

		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L, # intercept and treatment
				fit_fun = function(X_fit){
					start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit))
					if (estimate_only) {
						res = fast_logistic_regression_cpp(
							X_fit, private$y,
							start_beta = start_beta,
							smart_start = private$smart_default,
							optimization_alg = private$optimization_alg
						)
						list(b = res$b, ssq_b_2 = NA_real_)
					} else {
						fast_logistic_regression_with_var_cpp(
							X_fit, private$y,
							start_beta = start_beta,
							smart_start = private$smart_default,
							optimization_alg = private$optimization_alg
						)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
				}
			)

			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(attempt$fit$b, "beta")
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = match(2L, attempt$keep)
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
