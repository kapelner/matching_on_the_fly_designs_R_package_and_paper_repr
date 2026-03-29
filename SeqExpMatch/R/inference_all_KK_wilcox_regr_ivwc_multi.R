#' Multivariate Wilcox Rank-based Regression Compound Inference for KK Designs
#'
#' Fits a robust compound estimator for KK matching-on-the-fly designs using rank-based
#' regression (R-estimation) with the treatment indicator and all recorded covariates.
#' For matched pairs, it uses R-estimation on the within-pair differences of responses
#' and covariates. For reservoir subjects, it uses a standard rank-based linear model
#' including both treatment and covariates. This is the multivariate (covariate-adjusted)
#' variant.
#'
#' @export
InferenceAllKKWilcoxRegrMultiIVWC = R6::R6Class("InferenceAllKKWilcoxRegrMultiIVWC",
	lock_objects = FALSE,
	inherit = InferenceAllKKWilcoxRegrUnivIVWC,
	public = list(),
	private = list(
		include_covariates = function() TRUE
	)
)
