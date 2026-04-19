#' Robust-Regression IVWC Compound Inference for KK Designs
#'
#' Fits a variance-weighted compound estimator for KK matching-on-the-fly designs
#' with continuous responses using robust regression for matched-pair differences
#' and reservoir outcomes, with treatment and, optionally, all recorded covariates
#' as predictors.
#'
#' @export
InferenceContinKKRobustRegrIVWC = R6::R6Class("InferenceContinKKRobustRegrIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKRobustRegrIVWC,
	public = list(
	)
)
