#' Wald Incidence Inference
#'
#' Unadjusted incidence inference using the empirical risk difference together
#' with the standard unpooled Wald standard error and normal-approximation
#' confidence interval / hypothesis test.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidenceWald$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidenceWald = R6::R6Class("InferenceIncidenceWald",
	lock_objects = FALSE,
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		#' @description
		#' Initialize Wald incidence inference.
		#' @param des_obj A completed design object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @return A new \code{InferenceIncidenceWald} object.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
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
				self$compute_estimate()
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
