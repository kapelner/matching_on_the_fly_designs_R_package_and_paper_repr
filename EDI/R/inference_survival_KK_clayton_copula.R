#' Clayton Copula IVWC Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with survival responses
#' using a Clayton copula with Weibull AFT margins for matched pairs and standard
#' Weibull AFT regression for reservoir subjects.
#'
#' @export
InferenceIncidKKClaytonCopulaIVWC = R6::R6Class("InferenceIncidKKClaytonCopulaIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClaytonCopulaIVWC,
	public = list(
	)
)

#' Clayton Copula Combined-Likelihood Inference for KK Designs
#'
#' Fits a joint copula-based likelihood for KK matching-on-the-fly designs with
#' survival responses using a Clayton copula with Weibull AFT margins.
#'
#' @export
InferenceIncidKKClaytonCopulaOneLik = R6::R6Class("InferenceIncidKKClaytonCopulaOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClaytonCopulaOneLik,
	public = list(
	)
)
