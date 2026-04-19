#' Multivariate Robust-Regression Combined-Likelihood Inference for KK Designs
#'
#' Fits a stacked robust regression for KK matching-on-the-fly designs with
#' continuous responses using the treatment indicator and all recorded covariates.
#'
#' @export
InferenceContinMultiKKRobustRegrOneLik = R6::R6Class("InferenceContinMultiKKRobustRegrOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKRobustRegrOneLik,
	public = list(
	),
	private = list(
		include_covariates = function() TRUE
	)
)
