#' Multivariate Robust-Regression IVWC Compound Inference for KK Designs
#'
#' Fits a variance-weighted compound estimator for KK matching-on-the-fly designs
#' with continuous responses using robust regression for matched-pair differences
#' and reservoir outcomes, with treatment and all recorded covariates as predictors.
#'
#' @export
InferenceContinMultiKKRobustRegrIVWC = R6::R6Class("InferenceContinMultiKKRobustRegrIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKRobustRegrIVWC,
	public = list(
	),
	private = list(
		include_covariates = function() TRUE
	)
)
