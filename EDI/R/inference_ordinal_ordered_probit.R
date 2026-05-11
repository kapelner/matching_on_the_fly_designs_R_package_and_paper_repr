#' Ordered Probit Regression Inference for Ordinal Responses
#'
#' Fits an ordered probit regression for ordinal responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceOrdinalOrderedProbitRegr = R6::R6Class("InferenceOrdinalOrderedProbitRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize an ordered probit inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_default Whether to use smart optimizer start values by default.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
				assertFormula(model_formula, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_default = smart_default)
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
				X = matrix(private$w, ncol = 1L)
				colnames(X) = "treatment"
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(treatment = private$w, X_cov)
			}
			res = fast_ordinal_probit_regression_cpp(X = X, y = as.numeric(private$y))
			if (is.null(res) || length(res$b) < 1L || !is.finite(res$b[length(res$b)])){
				return(NA_real_)
			}
			as.numeric(res$b[length(res$b)])
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
			if (is.null(ctx)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			full_fit = list(params = ctx$full_params, neg_loglik = ctx$full_neg_loglik)
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = full_fit,
				fit_null = function(delta, start = NULL){
					res = tryCatch(
						fast_ordinal_probit_regression_cpp(
							X_fit, y,
							fixed_idx = j_treat, fixed_values = delta,
							start_params = start,
							smart_start = private$smart_default
						),
						error = function(e) NULL
					)
					if (is.null(res) || length(res) == 0) return(NULL)
					list(params = as.numeric(res$params), neg_loglik = as.numeric(res$neg_loglik))
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
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + length(sort(unique(private$y))) - 1L)
					if (estimate_only) {
						res = fast_ordinal_probit_regression_cpp(
							X_fit, private$y,
							start_params = start_params,
							smart_start = private$smart_default
						)
						if (is.null(res) || length(res) == 0) return(NULL)
						list(b = res$b, ssq_b_j = NA_real_, params = res$params, neg_loglik = res$neg_loglik)
					} else {
						res = fast_ordinal_probit_regression_with_var_cpp(
							X_fit, private$y,
							start_params = start_params,
							smart_start = private$smart_default
						)
						if (is.null(res) || length(res$b) == 0 || is.na(res$b[1])) return(NULL)
						list(b = res$b, ssq_b_j = res$ssq_b_j, params = res$params, neg_loglik = res$neg_loglik)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = length(mod$b)
					if (is.null(mod) || j_treat < 1L || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					ssq = mod$ssq_b_j
					!is.null(ssq) && is.finite(ssq) && ssq > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(attempt$fit$params, "params")
				private$best_X_colnames = setdiff(colnames(attempt$X), "treatment")
				n_alpha = length(attempt$fit$params) - ncol(attempt$X)
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = as.integer(n_alpha + 1L),
					full_params = attempt$fit$params,
					full_neg_loglik = attempt$fit$neg_loglik
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
