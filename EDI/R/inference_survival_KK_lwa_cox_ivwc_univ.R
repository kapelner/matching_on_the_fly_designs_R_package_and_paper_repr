#' Univariate LWA-style Marginal Cox Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with survival responses
#' using only the treatment indicator. For matched pairs, it uses a marginal Cox
#' proportional hazards model with Lee-Wei-Amato style cluster-robust variance. For
#' reservoir subjects, it uses standard Cox regression. Both models estimate log-hazard
#' ratios, which are then combined via a variance-weighted linear combination.
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
#' infer <- InferenceSurvivalUnivKKLWACoxIVWC$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceSurvivalUnivKKLWACoxIVWC = R6::R6Class("InferenceSurvivalUnivKKLWACoxIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKLWACoxIVWC,
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
		include_covariates = function() FALSE
	)
)
