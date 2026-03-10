#' Quantile Regression Combined-Likelihood Compound Estimator for KK Designs (stub, proportion)
#'
#' @description
#' Stub placeholder. Methods are not yet implemented.
#'
#' @export
SeqDesignInferencePropMultiKKQuantileRegrCombinedLikelihood = R6::R6Class("SeqDesignInferencePropMultiKKQuantileRegrCombinedLikelihood",
	inherit = SeqDesignInferenceAbstractKKQuantileRegrCombinedLikelihood,
	public = list(
		#' @description	Initialize the inference object.
		#' @param	seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param	tau				The quantile level on the logit scale, strictly between 0 and 1. Default is 0.5.
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, tau = 0.5, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "proportion")
			super$initialize(seq_des_obj, tau, qlogis, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		}
	)
)
