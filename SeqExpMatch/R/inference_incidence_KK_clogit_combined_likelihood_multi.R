#' Multivariate Conditional Logistic Combined-Likelihood Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with binary (incidence)
#' responses using the treatment indicator and all recorded covariates. Uses the combined
#' logistic likelihood over discordant matched-pair differences and reservoir subjects.
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
#'   response_type = "incidence",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 0, 1, 0, 1, 1, 0))
#' infer <- InferenceIncidMultiKKClogitCombinedLikelihood$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceIncidMultiKKClogitCombinedLikelihood = R6::R6Class("InferenceIncidMultiKKClogitCombinedLikelihood",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitCombinedLikelihood,
	public = list(



	),
	private = list(
		include_covariates = function() TRUE
	)
)
