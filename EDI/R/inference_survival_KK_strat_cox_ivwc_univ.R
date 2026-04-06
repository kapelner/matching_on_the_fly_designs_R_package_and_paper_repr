#' Univariate Stratified Cox / Standard Cox Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with survival responses
#' using only the treatment indicator. For matched pairs, it uses stratified Cox
#' proportional hazards regression (each pair is a stratum). For reservoir subjects, it
#' uses standard Cox regression. Both models estimate log-hazard ratios, which are then
#' combined via a variance-weighted linear combination.
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
#' infer <- InferenceSurvivalUnivKKStratCoxIVWC$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceSurvivalUnivKKStratCoxIVWC = R6::R6Class("InferenceSurvivalUnivKKStratCoxIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKStratCoxIVWC,
	public = list(




	),
	private = list(
		include_covariates = function() FALSE
	)
)
