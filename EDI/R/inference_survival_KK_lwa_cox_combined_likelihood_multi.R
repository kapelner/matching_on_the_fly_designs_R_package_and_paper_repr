#' Multivariate LWA-style Marginal Cox Combined-Likelihood Inference for KK Designs
#'
#' Fits a single joint marginal Cox model on all KK design data using the treatment
#' indicator and all recorded covariates, with Lee-Wei-Amato style cluster-robust
#' variance. Matched subjects share their pair ID as a cluster and reservoir subjects
#' are treated as singleton clusters.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneKK14$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "survival",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(
#'   c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
#'   c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' infer <- InferenceSurvivalMultiKKLWACoxCombinedLikelihood$
#'   new(seq_des, verbose =
#' FALSE)
#' infer
#'
InferenceSurvivalMultiKKLWACoxCombinedLikelihood = R6::R6Class("InferenceSurvivalMultiKKLWACoxCombinedLikelihood",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKLWACoxCombinedLikelihood,
	public = list(
		#' @description
		#' Initialize the Inference object.
		#'
		#' @param des_obj The design object.
		#' @param verbose If TRUE, print additional information.
		initialize = function(des_obj, verbose = FALSE) {
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj, verbose)
		}



	),
	private = list(
		include_covariates = function() TRUE
	)
)
