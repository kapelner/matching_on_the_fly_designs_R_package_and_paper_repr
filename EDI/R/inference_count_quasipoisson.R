#' Quasi-Poisson Regression Inference for Count Responses
#'
#' Fits a quasi-Poisson log-link regression for count responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountQuasiPoisson$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountQuasiPoisson = R6::R6Class("InferenceCountQuasiPoisson",
	lock_objects = FALSE,
	inherit = InferenceCountCompositeLikelihood,
	public = list(
				
		#' @description
		#' Initialize a quasi-Poisson regression inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},

		#' @description
		#' Computes the quasi-Poisson-regression estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description Indicates whether likelihood-based tests are supported.
		#' @return Always TRUE because companion Poisson tests are available.
		supports_likelihood_tests = function(){
			TRUE
		}
	),

	private = list(
		best_X_colnames = NULL,

		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			X
		},

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

			res = tryCatch(fast_poisson_regression_cpp(X = X, y = as.numeric(private$y)), error = function(e) NULL)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			as.numeric(res$b[2])
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		generate_mod = function(estimate_only = FALSE){
			# Use the common GLM fitting pattern
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = private$build_design_matrix(),
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					if (estimate_only) {
						res = fast_poisson_regression_cpp(X = X_fit, y = private$y)
						list(b = res$b, ssq_b_j = NA_real_, j_treat = j_treat)
					} else {
						res = fast_quasipoisson_regression_with_var_cpp(X = X_fit, y = private$y, j = j_treat)
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
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = which(attempt$keep == 2L)
				)
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		},

		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			companion_fit = tryCatch(
				fast_poisson_regression_with_var_cpp(
					X = X_fit,
					y = y,
					j = j_treat,
					start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
					smart_start = private$smart_default
				),
				error = function(e) NULL
			)
			if (is.null(companion_fit) || is.null(companion_fit$b) || length(companion_fit$b) < 2L || !is.finite(companion_fit$b[2])) return(NULL)
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = companion_fit,
				fit_null = function(delta, start = NULL){
					fast_poisson_regression_with_var_cpp(
						X = X_fit,
						y = y,
						j = j_treat,
						start_beta = start %||% private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
						fixed_idx = j_treat,
						fixed_values = delta,
						smart_start = private$smart_default
					)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					get_poisson_regression_score_cpp(X_fit, y, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				neg_loglik = function(fit){
					eta = as.numeric(X_fit %*% as.numeric(fit$b))
					-sum(y * eta - exp(eta) - lgamma(y + 1))
				}
			)
		}
	)
)
