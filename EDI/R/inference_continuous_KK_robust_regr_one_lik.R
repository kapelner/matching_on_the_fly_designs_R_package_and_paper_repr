#' Robust-Regression Combined-Likelihood Inference for KK Designs
#'
#' Fits a single stacked robust regression over matched-pair differences and
#' reservoir observations for KK matching-on-the-fly designs with continuous
#' responses, using the treatment indicator and, optionally, all recorded
#' covariates as predictors.
#'
#' @export
InferenceContinKKRobustRegrOneLik = R6::R6Class("InferenceContinKKRobustRegrOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKRobustRegrOneLik,
	public = list(
	)
)
