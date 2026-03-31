#' Simple Incidence Proportion Difference Inference
#'
#' Unadjusted incidence inference using the simple treated-minus-control
#' proportion difference with pooled equal-variance t inference.
#'
#' @export
InferenceIncidenceSimplePropDiffPooled = R6::R6Class("InferenceIncidenceSimplePropDiffPooled",
	lock_objects = FALSE,
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "incidence")
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
		}
	),

	private = list(
		get_standard_error = function(){
			if (is.null(private$cached_values$incidence_simple_prop_diff_se)) {
				private$compute_incidence_simple_prop_diff_pooled_components()
			}
			private$cached_values$incidence_simple_prop_diff_se
		},

		get_degrees_of_freedom = function(){
			if (is.null(private$cached_values$incidence_simple_prop_diff_df)) {
				private$compute_incidence_simple_prop_diff_pooled_components()
			}
			private$cached_values$incidence_simple_prop_diff_df
		},

		compute_incidence_simple_prop_diff_pooled_components = function(){
			if (is.null(private$cached_values$beta_hat_T)) {
				self$compute_treatment_estimate()
			}

			y_t = private$cached_values$yTs
			y_c = private$cached_values$yCs
			n_t = length(y_t)
			n_c = length(y_c)

			if (n_t <= 1L || n_c <= 1L) {
				private$cached_values$incidence_simple_prop_diff_se = NA_real_
				private$cached_values$incidence_simple_prop_diff_df = NA_real_
				return(invisible(NULL))
			}

			s2_t = stats::var(y_t)
			s2_c = stats::var(y_c)
			df = n_t + n_c - 2L
			s2_pooled = ((n_t - 1L) * s2_t + (n_c - 1L) * s2_c) / df
			var_hat = s2_pooled * (1 / n_t + 1 / n_c)

			private$cached_values$incidence_simple_prop_diff_se =
				if (is.finite(var_hat) && var_hat >= 0) sqrt(var_hat) else NA_real_

			private$cached_values$incidence_simple_prop_diff_df =
				if (is.finite(var_hat) && is.finite(df) && df > 0) {
					as.numeric(df)
				} else {
					NA_real_
				}
			invisible(NULL)
		}
	)
)
