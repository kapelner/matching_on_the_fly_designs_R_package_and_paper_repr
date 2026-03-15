#' Univariate Hurdle Negative Binomial Regression Inference for Count Responses
#'
#' @description
#' Fits a hurdle negative binomial regression for count responses using only the
#' treatment indicator in both the hurdle and count components. The hurdle
#' indicator model is fit on all subjects, and the zero-truncated negative
#' binomial count component is optimized in C++ on the positive-count subjects.
#' The reported treatment effect is the treatment coefficient from the conditional
#' count component, on the log-rate scale.
#'
#' @export
SeqDesignInferenceCountUnivHurdleNegBinRegr = R6::R6Class("SeqDesignInferenceCountUnivHurdleNegBinRegr",
	inherit = SeqDesignInferenceCountHurdleNegBinAbstract,
	private = list(
		hurdle_description = function() "Hurdle negative binomial regression"
	)
)
