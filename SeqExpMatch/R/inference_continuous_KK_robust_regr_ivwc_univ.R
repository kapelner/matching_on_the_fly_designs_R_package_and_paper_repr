#' Univariate Robust-Regression IVWC Compound Inference for KK Designs
#'
#' @description
#' Fits a variance-weighted compound estimator for KK matching-on-the-fly designs
#' with continuous responses using robust regression for matched-pair differences
#' and reservoir outcomes, with treatment as the only predictor.
#'
#' @export
InferenceContinUnivKKRobustRegrIVWC = R6::R6Class("InferenceContinUnivKKRobustRegrIVWC",
	inherit = InferenceAbstractKKRobustRegrIVWC,
	public = list(
	),
	private = list(
		include_covariates = function() FALSE
	)
)
