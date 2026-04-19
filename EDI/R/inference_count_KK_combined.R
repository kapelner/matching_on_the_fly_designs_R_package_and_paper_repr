#' KK Hurdle Poisson Combined-Likelihood Inference for Count Responses
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with count
#' responses using a joint likelihood over all subjects.
#'
#' @export
InferenceCountKKHurdlePoissonOneLik = R6::R6Class("InferenceCountKKHurdlePoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKHurdlePoissonOneLik,
	public = list(
	)
)

#' KK Poisson-Conditional-Poisson Combined-Likelihood Inference for Count Responses
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with count
#' responses using a joint likelihood over all subjects.
#'
#' @export
InferenceCountKKCPoissonOneLik = R6::R6Class("InferenceCountKKCPoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKPoissonCPoissonOneLik,
	public = list(
	)
)
