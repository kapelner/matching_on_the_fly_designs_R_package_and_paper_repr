#' Multivariate Robust-Regression Combined-Likelihood Inference for KK Designs
#'
#' Fits a stacked robust regression for KK matching-on-the-fly designs with
#' continuous responses using the treatment indicator and all recorded covariates.
#'
#' @export
InferenceContinMultiKKRobustRegrCombinedLikelihood = R6::R6Class("InferenceContinMultiKKRobustRegrCombinedLikelihood",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKRobustRegrCombinedLikelihood,
	public = list(
	),
	private = list(
		include_covariates = function() TRUE
	)
)
