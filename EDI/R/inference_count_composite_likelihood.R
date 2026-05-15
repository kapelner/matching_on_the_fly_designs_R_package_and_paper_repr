#' Count Composite Likelihood Inference Base
#'
#' Shared branch for count models that use a companion likelihood for testing
#' while the reported estimator may be robust or quasi-likelihood based.
#'
#' @keywords internal
InferenceCountCompositeLikelihood = R6::R6Class("InferenceCountCompositeLikelihood",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = list(
		#' @description Computes the treatment estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		generate_mod = function(estimate_only = FALSE) stop(class(self)[1], " must implement generate_mod()"),

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			
			model_output = private$generate_mod(estimate_only = estimate_only)
			private$cached_mod = model_output$mod %||% model_output
			
			if (is.null(model_output)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}
			
			private$cached_values$beta_hat_T = model_output$beta_hat_T %||% model_output$b[2]
			
			if (!is.null(model_output$b)) {
				private$set_fit_warm_start(
					as.numeric(model_output$b),
					type = "beta",
					fisher = model_output$XtWX %||% model_output$fisher_information
				)
			}

			if (estimate_only) return(invisible(NULL))
			
			ssq = model_output$ssq_b_j %||% model_output$ssq_b_2
			if (!is.null(ssq) && is.finite(ssq) && ssq > 0) {
				private$cached_values$s_beta_hat_T = sqrt(ssq)
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$df = model_output$df %||% Inf
		},

		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df %||% Inf
		},

		supports_likelihood_tests = function(){
			TRUE
		},

		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx)) return(NULL)
			
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			
			# Companion Poisson fit for testing
			companion_fit = tryCatch(
				fast_poisson_regression_with_var_cpp(
					X = X_fit,
					y = y,
					j = j_treat,
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
					smart_start = private$smart_default,
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit))
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
						warm_start_beta = start %||% private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
						warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit)),
						fixed_idx = j_treat,
						fixed_values = delta,
						smart_start = private$smart_default,
						optimization_alg = private$optimization_alg %||% "lbfgs"
					)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					as.numeric(fit$score %||% get_poisson_regression_score_cpp(X_fit, y, as.numeric(fit$b)))
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
