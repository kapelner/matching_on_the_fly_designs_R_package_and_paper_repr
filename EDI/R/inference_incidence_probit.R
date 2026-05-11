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
		#' @param smart_default Whether to use smart optimizer start values by default.
		#' @param optimization_alg  Optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_default = TRUE, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
				assertFormula(model_formula, null.ok = TRUE)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE, default = "newton_raphson")
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
				X = matrix(private$w, ncol = 1L)
				colnames(X) = "treatment"
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(treatment = private$w, X_cov)
			}
			start_params = private$get_fit_warm_start_for_length("params", ncol(X) + 1L)
			res = fast_ordinal_probit_regression_cpp(
				X = X,
				y = as.numeric(private$y),
				start_params = start_params,
				smart_start = private$smart_default,
				optimization_alg = private$optimization_alg
			)
			if (is.null(res) || length(res$b) < 1L || !is.finite(res$b[1L])){
				return(NA_real_)
			}
			private$set_fit_warm_start(as.numeric(res$params), "params")
			as.numeric(res$b[1L])
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
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = full_fit,
				fit_null = function(delta, start = NULL){
					res = tryCatch(
						fast_ordinal_probit_regression_cpp(
							X_fit,
							y,
							start_params = start,
							smart_start = private$smart_default,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						),
						error = function(e) NULL
					)
					if (is.null(res) || length(res) == 0L) return(NULL)
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
					start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + 1L)
					if (estimate_only) {
						res = fast_ordinal_probit_regression_cpp(
							X = X_fit,
							y = as.numeric(private$y),
							start_params = start_params,
							smart_start = private$smart_default,
							optimization_alg = private$optimization_alg
						)
						if (is.null(res) || length(res) == 0L) return(NULL)
						list(
							b = as.numeric(res$b),
							ssq_b_j = NA_real_,
							j_treat = j_treat,
							params = as.numeric(res$params),
							neg_loglik = as.numeric(res$neg_loglik)
						)
					} else {
						res = fast_ordinal_probit_regression_with_var_cpp(
							X = X_fit,
							y = as.numeric(private$y),
							start_params = start_params,
							smart_start = private$smart_default,
							optimization_alg = private$optimization_alg
						)
						if (is.null(res) || length(res$b) == 0L || is.na(res$b[1L])) return(NULL)
						list(
							b = as.numeric(res$b),
							ssq_b_j = as.numeric(res$ssq_b_j),
							j_treat = j_treat,
							params = as.numeric(res$params),
							neg_loglik = as.numeric(res$neg_loglik)
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
