#' Univariate Hurdle Poisson Regression Inference for Count Responses
#'
#' @description
#' Fits a hurdle Poisson regression for count responses using only the treatment
#' indicator in both the count and hurdle components. The reported treatment effect
#' is the treatment coefficient from the conditional truncated-Poisson count
#' component, on the log-rate scale.
#'
#' @export
SeqDesignInferenceCountUnivHurdlePoissonRegr = R6::R6Class("SeqDesignInferenceCountUnivHurdlePoissonRegr",
	inherit = SeqDesignInferenceCountZeroAugmentedPoissonAbstract,
	private = list(
		za_family = function() glmmTMB::truncated_poisson(link = "log"),
		za_description = function() "Hurdle Poisson regression"
	)
)
