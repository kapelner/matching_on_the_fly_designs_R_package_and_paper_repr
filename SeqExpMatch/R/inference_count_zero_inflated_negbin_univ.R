#' Univariate Zero-Inflated Negative Binomial Regression Inference for Count Responses
#'
#' @description
#' Fits a zero-inflated negative binomial regression for count responses using only
#' the treatment indicator in both the count and zero-inflation components. The
#' reported treatment effect is the treatment coefficient from the conditional
#' negative binomial count component, on the log-rate scale.
#'
#' @export
SeqDesignInferenceCountUnivZeroInflatedNegBinRegr = R6::R6Class("SeqDesignInferenceCountUnivZeroInflatedNegBinRegr",
	inherit = SeqDesignInferenceCountZeroAugmentedPoissonAbstract,
	private = list(
		za_family = function() glmmTMB::nbinom2(link = "log"),
		za_description = function() "Zero-inflated negative binomial regression"
	)
)
