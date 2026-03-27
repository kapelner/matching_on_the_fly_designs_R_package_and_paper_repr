#' Univariate Wilcox Rank-based Regression Compound Inference for KK Designs
#'
#' @description
#' Fits a robust compound estimator for KK matching-on-the-fly designs using rank-based
#' regression (R-estimation) with only the treatment indicator (no additional covariates).
#' For matched pairs, it uses R-estimation on the within-pair response differences.
#' For reservoir subjects, it uses a standard rank-based linear model with the treatment
#' assignment as the sole predictor. This is the univariate (covariate-free) variant.
#'
#' @export
InferenceAllKKWilcoxRegrUnivIVWC = R6::R6Class("InferenceAllKKWilcoxRegrUnivIVWC",
	inherit = InferenceAbstractKKWilcoxRegrIVWC,
	public = list(

	),
	private = list(
		include_covariates = function() FALSE
	)
)
