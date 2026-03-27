#' Exact Zhang Incidence Inference for KK or Bernoulli Designs
#'
#' @description
#' Implements the exact test-inversion inference of Zhang (2026). This method
#' provides a combined confidence interval and p-value by partitioning the
#' data into matched pairs (analyzed via a binomial test) and a reservoir of
#' unmatched subjects (analyzed via Fisher's Exact Test).
#'
#' @details
#' This method is specifically designed for binary (incidence) outcomes in
#' KK matching-on-the-fly designs or Bernoulli designs. It uses a bisection
#' algorithm to invert the combined p-value for confidence interval construction.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceExact methods
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
DesignInferenceIncidExactZhang = R6::R6Class("DesignInferenceIncidExactZhang",
	inherit = DesignInferenceExact,
	public = list(
		#' @description
		#' Initialize the Zhang inference object.
		#' @param des_obj A completed \code{Design} object (incidence response).
		#' @param num_cores Number of CPU cores for bisection search.
		#' @param verbose Flag for progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "incidence")
			super$initialize(des_obj, num_cores, verbose)
		},

		#' @description
		#' Returns the incidence treatment estimate (log odds ratio).
		compute_treatment_estimate = function(){
			stats = private$get_exact_zhang_stats()
			zhang_incid_treatment_estimate(stats)
		}
	)
)
