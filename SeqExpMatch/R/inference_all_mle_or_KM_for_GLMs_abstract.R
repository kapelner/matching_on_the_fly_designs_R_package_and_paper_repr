#' GLM and Kaplan-Meier Inference
#'
#' @description
#' Abstract class providing MLE/KM-based inference methods for GLM and survival models.
#'
#' @keywords internal
InferenceMLEorKMforGLMs = R6::R6Class("InferenceMLEorKMforGLMs",
	inherit = InferenceAsymp,
	public = list(

		# @description
		# Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		# @param des_obj		A DesignSeqOneByOne object whose entire n subjects are assigned and response y is recorded within.
		# @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		# 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference
		# 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs),
		# 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		# @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}

		#' @description
		#' Computes the treatment estimate using the underlying model.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		}
	),

	private = list(
		generate_mod = function() stop(class(self)[1], " must implement generate_mod()"),

		get_standard_error = function(){
			private$shared()
			private$cached_values$s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			private$shared()
			private$cached_values$df
		},

		shared = function(){
			if (!is.null(private$cached_values$is_z)) return(invisible(NULL))
			model_output = private$generate_mod() #abstract function implemented by daughter classes. Should return a list with 'b' and 'ssq_b_2'.
			private$cached_values$beta_hat_T = model_output$b[2]

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
