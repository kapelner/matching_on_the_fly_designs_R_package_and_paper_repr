#' Univariate Conditional Logistic Combined-Likelihood Compound Inference for KK Designs
#'
#' @description
#' Fits a compound estimator for KK matching-on-the-fly designs with binary (incidence)
#' responses using only the treatment indicator (no additional covariates). Uses the
#' combined logistic likelihood over discordant matched-pair differences and reservoir
#' subjects.
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
#' infer <- InferenceIncidUnivKKClogitCombinedLikelihood$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceIncidUnivKKClogitCombinedLikelihood = R6::R6Class("InferenceIncidUnivKKClogitCombinedLikelihood",
	inherit = InferenceAbstractKKClogitCombinedLikelihood,
	public = list(



	),
	private = list(
		include_covariates = function() FALSE
	)
)
