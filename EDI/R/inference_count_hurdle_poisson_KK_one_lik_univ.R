#' Univariate KK Hurdle Poisson Combined-Likelihood Inference for Count Responses
#'
#' Fits a KK hurdle-Poisson combined-likelihood model for count responses using
#' only the treatment indicator in both the hurdle and positive-count
#' components. Matched pairs contribute pair-specific random intercepts, while
#' reservoir subjects contribute through a common fixed reservoir intercept. The
#' reported treatment effect is the treatment coefficient from the positive-count
#' component on the log-rate scale.
#'
#' @export
InferenceCountUnivKKHurdlePoissonOneLik = R6::R6Class("InferenceCountUnivKKHurdlePoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKHurdlePoissonOneLik,
	private = list(
		include_covariates = function() FALSE
	)
)
