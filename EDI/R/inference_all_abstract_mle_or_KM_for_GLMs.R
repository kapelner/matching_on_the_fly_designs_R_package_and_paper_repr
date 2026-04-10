#' GLM and Kaplan-Meier Inference
#'
#' Abstract class providing MLE/KM-based inference methods for GLM and survival models.
#'
#' @keywords internal
InferenceMLEorKMforGLMs = R6::R6Class("InferenceMLEorKMforGLMs",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Computes the treatment estimate using the underlying model.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		}
	),

	private = list(
		generate_mod = function(estimate_only = FALSE) stop(class(self)[1], " must implement generate_mod()"),

		create_bootstrap_worker_state = function(){
			private$create_design_backed_bootstrap_worker_state()
		},

		load_bootstrap_sample_into_worker = function(worker_state, indices){
			private$load_bootstrap_sample_into_design_backed_worker(worker_state, indices)
		},

		compute_bootstrap_worker_estimate = function(worker_state){
			private$compute_bootstrap_worker_estimate_via_compute_treatment_estimate(worker_state)
		},

		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			has_cached_se = !is.null(private$cached_values$s_beta_hat_T) &&
				length(private$cached_values$s_beta_hat_T) == 1L &&
				isTRUE(is.finite(private$cached_values$s_beta_hat_T))
			if (isTRUE(!is.null(private$cached_values$beta_hat_T) && (estimate_only || has_cached_se))) return(invisible(NULL))
			model_output = private$generate_mod(estimate_only = estimate_only) #abstract function implemented by daughter classes. Should return a list with 'b' and 'ssq_b_2'.
			private$cached_mod = model_output
			private$cached_values$beta_hat_T = model_output$b[2]
			if (estimate_only) return(invisible(NULL))

			ssq = model_output$ssq_b_2
			if (!is.null(ssq) && !is.na(ssq) && ssq > 0) {
				private$cached_values$s_beta_hat_T = sqrt(ssq)
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$is_z = TRUE # This remains true for asymptotic inference
			private$cached_values$df = model_output$df %||% NA_real_
		}
	)
)
