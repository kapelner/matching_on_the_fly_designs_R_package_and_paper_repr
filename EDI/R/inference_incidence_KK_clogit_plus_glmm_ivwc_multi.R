#' Multivariate Conditional Logistic Plus GLMM IVWC Inference for KK Designs
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs. The treatment effect is estimated by the conditional-logistic
#' component; the GLMM component includes only the intercept and covariates.
#'
#' @export
InferenceIncidMultiKKClogitPlusGLMMIVWC = R6::R6Class("InferenceIncidMultiKKClogitPlusGLMMIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitPlusGLMM,
	public = list(),
	private = list(
		combine_reservoir_into_glmm = function() FALSE
	)
)
