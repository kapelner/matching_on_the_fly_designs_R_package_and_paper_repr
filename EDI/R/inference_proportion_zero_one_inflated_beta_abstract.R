#' Zero/One-Inflated Beta Inference for Proportion Responses
#'
#' Internal base class for non-KK zero/one-inflated beta regression models. The
#' response is modeled as a three-component mixture with point masses at 0 and 1
#' plus a beta-distributed interior component on \eqn{(0, 1)}. The reported
#' treatment effect is the treatment coefficient from the beta mean submodel, on
#' the logit scale.
#'
#' @details
#' The inflation probabilities for exact 0 and exact 1 are intercept-only in this
#' first implementation. The beta mean submodel uses treatment alone in the
#' univariate class and treatment plus covariates in the multivariate class.
#'
#' @keywords internal
#' @noRd
InferencePropZeroOneInflatedBetaAbstract = R6::R6Class("InferencePropZeroOneInflatedBetaAbstract",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a zero-one-inflated beta regression inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "proportion")
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
			X = if (length(X_cols) == 0L) {
				cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				cbind(`(Intercept)` = 1, treatment = private$w, X_data[, X_cols, drop = FALSE])
			}
			
			res = tryCatch(fast_zero_one_inflated_beta_cpp(X, private$y, init = rep(0, ncol(X) + 3)), error = function(e) NULL)
			if (is.null(res) || length(res$b) < 2 || !is.finite(res$b[2])) return(NA_real_)
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
			full_params = ctx$full_params
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta){
					init = if (!is.null(full_params)) as.numeric(full_params) else rep(0, ncol(X_fit) + 3)
					res = tryCatch(fast_zero_one_inflated_beta_cpp(X_fit, y, init = init,
					               fixed_idx = j_treat, fixed_values = delta), error = function(e) NULL)
					if (is.null(res)) return(NULL)
					res
				},
				score = function(fit){
					get_zero_one_inflated_beta_score_cpp(X_fit, y, as.numeric(fit$params))
				},
				observed_information = function(fit){
					-get_zero_one_inflated_beta_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				fisher_information = function(fit){
					-get_zero_one_inflated_beta_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				information = function(fit){
					-get_zero_one_inflated_beta_hessian_cpp(X_fit, y, as.numeric(fit$params))
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
					init = rep(0, ncol(X_fit) + 3)
					res = tryCatch(fast_zero_one_inflated_beta_cpp(X_fit, private$y, init = init),
					               error = function(e) NULL)
					if (is.null(res)) return(NULL)
					ssq_b_j = if (!is.null(res$vcov) && nrow(res$vcov) >= j_treat) {
						as.numeric(res$vcov[j_treat, j_treat])
					} else {
						NA_real_
					}
					list(b = as.numeric(res$b), ssq_b_j = ssq_b_j, j_treat = j_treat,
					     params = as.numeric(res$params), neg_loglik = as.numeric(res$neg_loglik))
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
					j_treat = attempt$fit$j_treat,
					full_params = attempt$fit$params
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		}
	)
)
