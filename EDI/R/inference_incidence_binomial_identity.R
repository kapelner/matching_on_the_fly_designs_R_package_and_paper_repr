#' Binomial Identity Risk Difference Inference for Incidence Responses
#'
#' Fits a binomial identity-link regression for binary (incidence) responses
#' using the treatment indicator and, optionally, all recorded covariates as
#' predictors. The treatment effect is reported as a risk difference.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidBinomialIdentityRiskDiff$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidBinomialIdentityRiskDiff = R6::R6Class("InferenceIncidBinomialIdentityRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
				
		#' @description Initialize a binomial identity-link risk-difference inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
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
				X = cbind(1, private$w)
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(1, treatment = private$w, X_cov)
			}
			res = tryCatch(fast_identity_binomial_regression_cpp(X = X, y = as.numeric(private$y)), error = function(e) NULL)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
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
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta){
					res = tryCatch(fast_identity_binomial_regression_cpp(X = X_fit, y = y,
					               fixed_idx = j_treat, fixed_values = delta), error = function(e) NULL)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					res
				},
				score = function(fit){
					get_identity_binomial_regression_score_cpp(X_fit, y, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_identity_binomial_regression_hessian_cpp(X_fit, y, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_identity_binomial_regression_hessian_cpp(X_fit, y, as.numeric(fit$b))
				},
				information = function(fit){
					-get_identity_binomial_regression_hessian_cpp(X_fit, y, as.numeric(fit$b))
				},
				neg_loglik = function(fit){
					mu = as.numeric(X_fit %*% as.numeric(fit$b))
					-sum(y * log(pmax(mu, 1e-15)) + (1 - y) * log(pmax(1 - mu, 1e-15)))
				}
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
					if (estimate_only) {
						res = tryCatch(fast_identity_binomial_regression_cpp(X = X_fit, y = private$y),
						               error = function(e) NULL)
						if (is.null(res)) return(NULL)
						list(b = res$b, ssq_b_j = NA_real_, j_treat = j_treat)
					} else {
						res = tryCatch(fast_identity_binomial_regression_with_var_cpp(X = X_fit, y = private$y, j = j_treat),						               error = function(e) NULL)
						if (is.null(res)) return(NULL)
						res$j_treat = j_treat
						res
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
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		}
	)
)
