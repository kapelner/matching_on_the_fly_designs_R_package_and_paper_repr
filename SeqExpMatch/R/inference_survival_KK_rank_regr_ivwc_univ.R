#' Univariate Survival Rank-based Regression (AFT) Compound Inference for KK Designs
#'
#' @description
#' Fits a robust compound estimator for KK matching-on-the-fly designs with survival
#' responses using rank-based estimating equations. For matched pairs, it uses a
#' rank-based AFT model with clustering. For reservoir subjects, it uses a standard
#' rank-based AFT model. Both models estimate log-time ratios, which are then combined
#' via a variance-weighted linear combination. This is the univariate (covariate-free)
#' variant.
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
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(
#'   c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
#'   c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' infer <- InferenceSurvivalUnivKKRankRegrIVWC$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceSurvivalUnivKKRankRegrIVWC = R6::R6Class("InferenceSurvivalUnivKKRankRegrIVWC",
	inherit = InferenceAbstractKKSurvivalRankRegrIVWC,
	public = list(




	),
	private = list(
		include_covariates = function() FALSE
	)
)
