#' Univariate Zero-Inflated Poisson Regression Inference for Count Responses
#'
#' @description
#' Fits a zero-inflated Poisson regression for count responses using only the
#' treatment indicator in both the count and zero-inflation components. The
#' reported treatment effect is the treatment coefficient from the conditional
#' Poisson count component, on the log-rate scale.
#'
#' @export
SeqDesignInferenceCountUnivZeroInflatedPoissonRegr = R6::R6Class("SeqDesignInferenceCountUnivZeroInflatedPoissonRegr",
	inherit = SeqDesignInferenceCountZeroAugmentedPoissonAbstract,
	private = list(
		za_family = function() poisson(link = "log"),
		za_description = function() "Zero-inflated Poisson regression"
	)
)
