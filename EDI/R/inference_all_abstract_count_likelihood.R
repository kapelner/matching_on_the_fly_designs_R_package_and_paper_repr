#' Count-Specific Likelihood Inference
#'
#' @name InferenceCountLikelihood
#' @description Intermediate base class for count-based likelihood families
#' (Poisson, Negative Binomial, Zero-Inflated, Hurdle). This class centralizes
#' count-specific parameter packing, warm starts, and likelihood dispatch.
#'
#' @keywords internal
InferenceCountLikelihood = R6::R6Class("InferenceCountLikelihood",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = list(
		#' @description Computes the treatment estimate using the underlying model.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		# --- Count-specific shared logic ---

		# 1. Parameter packing and unpacking
		# Count models often pack coefficients (b) and auxiliary parameters (e.g. log-theta).
		# We use "params" as the standard tag for packed parameters in count models.

		# 2. Warm starts
		# Centralizes the use of "params" and "beta" tags for count warm starts.

		# 3. Estimate, SE, and df caching
		# Re-implements the basic caching pattern from StdModCache but scoped to counts.
		
		generate_mod = function(estimate_only = FALSE) stop(class(self)[1], " must implement generate_mod()"),

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			
			model_output = private$generate_mod(estimate_only = estimate_only)
			private$cached_mod = model_output
			
			if (is.null(model_output)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}
			
			# Count models usually return treatment effect at index 2 of the cond component
			# or in a specific field. We allow models to set beta_hat_T directly in output.
			private$cached_values$beta_hat_T = model_output$beta_hat_T %||% model_output$b[2]
			
			if (!is.null(model_output$params) || !is.null(model_output$b)) {
				private$set_fit_warm_start(
					as.numeric(model_output$params %||% model_output$b),
					type = if (!is.null(model_output$params)) "params" else "beta",
					fisher = model_output$fisher_information %||% model_output$XtWX
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
			# Try information-based SE first if supported
			if (isTRUE(private$supports_information_preference())) {
				se = tryCatch(private$compute_standard_error_from_information_matrix(), error = function(e) NA_real_)
				if (is.finite(se)) return(se)
			}
			private$cached_values$s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df %||% Inf
		},
		get_backend_warm_start_args = function(expected_length, expected_fisher_dim = expected_length) {
			private$get_optimal_warm_start_config(expected_length, expected_fisher_dim)
		},

		# --- Likelihood test support ---

		supports_likelihood_tests = function(){
			TRUE
		},

		get_likelihood_test_spec = function(){
			# This is still abstract, but we provide the structure
			NULL
		},

		compute_score_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "score")
		},
		compute_gradient_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "gradient")
		},
		compute_lik_ratio_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio")
		},
		compute_score_confidence_interval_impl = function(alpha){
			private$invert_test_pval_confidence_interval(alpha)
		},
		compute_gradient_confidence_interval_impl = function(alpha){
			private$invert_test_pval_confidence_interval(alpha)
		},
		compute_lik_ratio_confidence_interval_impl = function(alpha){
			private$invert_test_pval_confidence_interval(alpha)
		}
	)
)
