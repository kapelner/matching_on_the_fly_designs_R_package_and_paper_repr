#' Poisson Regression Inference for Count Responses
#'
#' Fits a Poisson log-link regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceCountPoisson = R6::R6Class("InferenceCountPoisson",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(
				
		#' @description
		#' Initialize a Poisson regression inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose			Whether to print progress messages.
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
		#' Computes the Poisson-regression estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		}
	),

	private = list(
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}

			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			
			if (length(Xmm_cols) == 0L){
				Xmm = cbind(1, private$w)
			} else {
				X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				Xmm = cbind(1, treatment = private$w, X_cov)
			}

			res = tryCatch(fast_poisson_regression_cpp(X = Xmm, y = as.numeric(private$y)), error = function(e) NULL)
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
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta){
					fast_poisson_regression_with_var_cpp(
						Xmm = X_fit,
						y = y,
						j = j_treat,
						fixed_idx = j_treat,
						fixed_values = delta
					)
				},
				score = function(fit){
					get_poisson_regression_score_cpp(X_fit, y, as.numeric(fit$b))
				},
				information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				neg_loglik = function(fit){
					eta = as.numeric(X_fit %*% as.numeric(fit$b))
					-sum(y * eta - exp(eta) - lgamma(y + 1))
				}
			)
		},

		generate_mod = function(estimate_only = FALSE){
			# Use the common GLM fitting pattern
			X_data = private$get_X()
			if (is.null(X_data) || ncol(X_data) == 0) {
				X_full = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X_full = cbind(`(Intercept)` = 1, treatment = private$w, X_data)
			}
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					if (estimate_only) {
						res = fast_poisson_regression_cpp(X = X_fit, y = private$y)
						list(b = res$b, ssq_b_j = NA_real_, j_treat = j_treat)
					} else {
						res = fast_poisson_regression_with_var_cpp(Xmm = X_fit, y = private$y, j = j_treat)
						res$j_treat = j_treat
						res
					}
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod)) return(FALSE)
					j_treat = mod$j_treat
					if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)

			if (!is.null(attempt$fit)){
				private$best_Xmm_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = which(attempt$keep == 2L)
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			w_mat = permutations$w_mat
			if (is.null(w_mat)) return(NULL)
			X_covars = private$X
			log_transform = transform_responses == "log"
			compute_poisson_distr_parallel_cpp(as.numeric(y), X_covars, w_mat, as.numeric(delta), log_transform, private$n_cpp_threads(ncol(w_mat)))
		}
	)
)
