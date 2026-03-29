#' Univariate Robust-Regression Combined-Likelihood Inference for KK Designs
#'
#' Fits a stacked robust regression for KK matching-on-the-fly designs with
#' continuous responses using only the treatment indicator.
#'
#' @export
InferenceContinUnivKKRobustRegrCombinedLikelihood = R6::R6Class("InferenceContinUnivKKRobustRegrCombinedLikelihood",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKRobustRegrCombinedLikelihood,
	public = list(
	),
	private = list(
		include_covariates = function() FALSE
	)
)
