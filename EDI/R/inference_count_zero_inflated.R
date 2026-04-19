#' Zero-Inflated Poisson Regression Inference for Count Responses
#'
#' Fits a zero-inflated Poisson regression for count responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceCountZeroInflatedPoisson = R6::R6Class("InferenceCountZeroInflatedPoisson",
	lock_objects = FALSE,
	inherit = InferenceCountZeroAugmentedPoissonAbstract,
	public = list(
	),
	private = list(
		za_family = function() stats::poisson(link = "log"),
		za_description = function() "Zero-Inflated Poisson"
	)
)

#' Zero-Inflated Negative Binomial Regression Inference for Count Responses
#'
#' Fits a zero-inflated negative binomial regression for count responses using
#' the treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceCountZeroInflatedNegBin = R6::R6Class("InferenceCountZeroInflatedNegBin",
	lock_objects = FALSE,
	inherit = InferenceCountZeroAugmentedPoissonAbstract,
	public = list(
	),
	private = list(
		za_family = function() glmmTMB::nbinom2(link = "log"),
		za_description = function() "Zero-Inflated Negative Binomial"
	)
)
