#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#'
#'
#' @export
InferenceIncidUnivLogRegr = R6::R6Class("InferenceIncidUnivLogRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#'   The default is 1 for serial computation. For simple
		#'   estimators (e.g. mean difference and KK compound),
		#'   parallelization is achieved with zero-overhead C++ OpenMP.
		#'   For complex models (e.g. GLMs),
		#'   parallelization falls back to R's
		#'   \code{parallel::mclapply}, which incurs
		#'   session-forking overhead.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{TRUE}.
		#' @param make_fork_cluster Whether to use a fork cluster for parallelization.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "incidence")
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
		}
	),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			# Use the standard design matrix order: Intercept, Treatment, Covariates
			# This ensures treatment is at index 2 for fast_logistic_regression_with_var_cpp
			Xmm = cbind(1, private$w)
			colnames(Xmm) = c("(Intercept)", "treatment")

			if (estimate_only) {
				res = fast_logistic_regression_cpp(Xmm, private$y)
				list(b = res$b, ssq_b_2 = NA_real_)
			} else {
				# Call C++ function that computes both coefficients and variance
				fast_logistic_regression_with_var_cpp(Xmm, private$y)
			}
		}
	)
)
