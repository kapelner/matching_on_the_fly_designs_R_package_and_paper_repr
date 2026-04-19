#' Univariate Robust-Regression Combined-Likelihood Inference for KK Designs
#'
#' Fits a stacked robust regression for KK matching-on-the-fly designs with
#' continuous responses using only the treatment indicator.
#'
#' @export
InferenceContinUnivKKRobustRegrOneLik = R6::R6Class("InferenceContinUnivKKRobustRegrOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKRobustRegrOneLik,
	public = list(
	),
	private = list(
		include_covariates = function() FALSE
	)
)
