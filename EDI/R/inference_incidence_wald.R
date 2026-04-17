#' Wald Incidence Inference
#'
#' Unadjusted incidence inference using the empirical risk difference together
#' with the standard unpooled Wald standard error and normal-approximation
#' confidence interval / hypothesis test.
#'
#' @export
InferenceIncidenceWald = R6::R6Class("InferenceIncidenceWald",
	lock_objects = FALSE,
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		#' @description
		#' Initialize Wald incidence inference.
		#' @param des_obj A completed design object.
		#' @param verbose Whether to print progress messages.
		#' @return A new \code{InferenceIncidenceWald} object.
		initialize = function(des_obj, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$cached_values$is_z = TRUE
		}
	),

	private = list(
		get_standard_error = function(){
			if (is.null(private$cached_values$incidence_wald_se)) {
				private$compute_incidence_wald_components()
			}
			private$cached_values$incidence_wald_se
		},

		get_degrees_of_freedom = function(){
			NA_real_
		},

		compute_incidence_wald_components = function(){
			if (is.null(private$cached_values$beta_hat_T)) {
				self$compute_treatment_estimate()
			}

			y_t = private$cached_values$yTs
			y_c = private$cached_values$yCs
			n_t = length(y_t)
			n_c = length(y_c)

			if (n_t == 0L || n_c == 0L) {
				private$cached_values$incidence_wald_se = NA_real_
				return(invisible(NULL))
			}

			p_t = mean(y_t)
			p_c = mean(y_c)
			var_hat = p_t * (1 - p_t) / n_t + p_c * (1 - p_c) / n_c
			private$cached_values$incidence_wald_se =
				if (is.finite(var_hat) && var_hat >= 0) sqrt(var_hat) else NA_real_
			invisible(NULL)
		}
	)
)
