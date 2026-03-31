#' Quantile Regression Combined-Likelihood Compound Estimator for KK Designs (Continuous)
#'
#' Fits the combined stacked quantile regression (matched-pair differences + reservoir)
#' using the treatment indicator and all recorded covariates for continuous responses.
#' Minimises the joint check-function loss over both data sources simultaneously.
#' Inference is based on the stacked combined-likelihood quantile-regression fit.
#'
#' @export
InferenceContinMultKKQuantileRegrCombinedLikelihood = R6::R6Class("InferenceContinMultKKQuantileRegrCombinedLikelihood",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKQuantileRegrCombinedLikelihood,
	public = list(
		#' @description	Initialize the inference object.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param	tau				The quantile level for regression, strictly between 0 and 1. Default is 0.5.
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		#' @examples
		#' set.seed(1)
		#' x_dat <- data.frame(
		#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
		#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
		#' )
		#' seq_des <- DesignSeqOneByOneKK14$new(n = nrow(x_dat), response_type = "continuous", verbose =
		#' FALSE)
		#' for (i in seq_len(nrow(x_dat))) {
		#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
		#' }
		#' seq_des$add_all_subject_responses(c(1.2, 0.9, 1.5, 1.8, 2.1, 1.7, 2.6, 2.2))
		#' infer <- InferenceContinMultKKQuantileRegrCombinedLikelihood$new(seq_des, verbose
		#' = FALSE)
		#' infer
		#'
		initialize = function(des_obj, tau = 0.5, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "continuous")
			super$initialize(des_obj, tau, identity, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},


		#' @description
		#' Returns the estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE) super$compute_treatment_estimate()


	)
)
