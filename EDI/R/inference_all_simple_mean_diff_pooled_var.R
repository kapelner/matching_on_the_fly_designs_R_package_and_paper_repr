#' Simple Mean Difference Inference with Pooled Variance
#'
#' Unadjusted mean-difference inference using the simple treated-minus-control
#' difference with pooled equal-variance t inference.
#'
#' @export
InferenceAllSimpleMeanDiffPooledVar = R6::R6Class("InferenceAllSimpleMeanDiffPooledVar",
	lock_objects = FALSE,
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		#' @description
		#' Initialize simple pooled-variance inference.
		#' @param des_obj A completed design object.
		#' @param verbose Whether to print progress messages.
		#' @return A new \code{InferenceAllSimpleMeanDiffPooledVar} object.
		initialize = function(des_obj,  verbose = FALSE){
			super$initialize(des_obj, verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	),

	private = list(
		get_standard_error = function(){
			if (is.null(private$cached_values$simple_mean_diff_pooled_se)) {
				private$compute_simple_mean_diff_pooled_components()
			}
			private$cached_values$simple_mean_diff_pooled_se
		},

		get_degrees_of_freedom = function(){
			if (is.null(private$cached_values$simple_mean_diff_pooled_df)) {
				private$compute_simple_mean_diff_pooled_components()
			}
			private$cached_values$simple_mean_diff_pooled_df
		},

		compute_simple_mean_diff_pooled_components = function(){
			if (is.null(private$cached_values$beta_hat_T)) {
				self$compute_treatment_estimate()
			}

			y_t = private$cached_values$yTs
			y_c = private$cached_values$yCs
			n_t = length(y_t)
			n_c = length(y_c)

			if (n_t <= 1L || n_c <= 1L) {
				private$cached_values$simple_mean_diff_pooled_se = NA_real_
				private$cached_values$simple_mean_diff_pooled_df = NA_real_
				return(invisible(NULL))
			}

			s2_t = stats::var(y_t)
			s2_c = stats::var(y_c)
			df = n_t + n_c - 2L
			s2_pooled = ((n_t - 1L) * s2_t + (n_c - 1L) * s2_c) / df
			var_hat = s2_pooled * (1 / n_t + 1 / n_c)

			private$cached_values$simple_mean_diff_pooled_se =
				if (is.finite(var_hat) && var_hat >= 0) sqrt(var_hat) else NA_real_

			private$cached_values$simple_mean_diff_pooled_df =
				if (is.finite(var_hat) && is.finite(df) && df > 0) {
					as.numeric(df)
				} else {
					NA_real_
				}
			invisible(NULL)
		}
	)
)
