#' Hurdle Poisson Regression Inference for Count Responses
#'
#' Fits a hurdle Poisson regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceCountHurdlePoisson = R6::R6Class("InferenceCountHurdlePoisson",
	lock_objects = FALSE,
	inherit = InferenceCountZeroAugmentedPoissonAbstract,
	public = list(
	),
	private = list(
		za_family = function() glmmTMB::truncated_poisson(link = "log"),
		za_description = function() "Hurdle Poisson"
	)
)

#' Hurdle Negative Binomial Regression Inference for Count Responses
#'
#' Fits a hurdle negative binomial regression for count responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceCountHurdleNegBin = R6::R6Class("InferenceCountHurdleNegBin",
	lock_objects = FALSE,
	inherit = InferenceCountHurdleNegBinAbstract,
	public = list(
	),
	private = list(
		hurdle_description = function() "Hurdle Negative Binomial"
	)
)
