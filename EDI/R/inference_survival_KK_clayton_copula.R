#' Clayton Copula IVWC Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with survival responses
#' using a Clayton copula with Weibull AFT margins for matched pairs and standard
#' Weibull AFT regression for reservoir subjects.
#'
#' @examples
#' \dontrun{
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalKKClaytonCopulaOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#' @export
InferenceSurvivalKKClaytonCopulaIVWC = R6::R6Class("InferenceSurvivalKKClaytonCopulaIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClaytonCopulaIVWC,
	public = list(
	)
)
InferenceIncidKKClaytonCopulaIVWC = InferenceSurvivalKKClaytonCopulaIVWC

#' Clayton Copula Combined-Likelihood Inference for KK Designs
#'
#' Fits a joint copula-based likelihood for KK matching-on-the-fly designs with
#' survival responses using a Clayton copula with Weibull AFT margins.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalKKClaytonCopulaOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalKKClaytonCopulaOneLik = R6::R6Class("InferenceSurvivalKKClaytonCopulaOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClaytonCopulaOneLik,
	public = list(
	)
)
InferenceIncidKKClaytonCopulaOneLik = InferenceSurvivalKKClaytonCopulaOneLik
