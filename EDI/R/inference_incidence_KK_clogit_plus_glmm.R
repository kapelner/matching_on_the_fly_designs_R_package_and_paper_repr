#' Conditional Logistic Plus GLMM IVWC Inference for KK Designs
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs. The treatment effect is estimated by the conditional-logistic
#' component; the GLMM component includes only the intercept and covariates.
#'
#' @export
InferenceIncidKKClogitPlusGLMMIVWC = R6::R6Class("InferenceIncidKKClogitPlusGLMMIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitPlusGLMM,
	public = list(),
	private = list(
		combine_reservoir_into_glmm = function() FALSE
	)
)

#' Conditional Logistic Plus GLMM Combined-Likelihood Inference for KK Designs
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs. The treatment effect is estimated by the conditional-logistic
#' component; the GLMM component includes only the intercept and covariates.
#'
#' @export
InferenceIncidKKClogitPlusGLMMOneLik = R6::R6Class("InferenceIncidKKClogitPlusGLMMOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitPlusGLMM,
	public = list(),
	private = list(
		combine_reservoir_into_glmm = function() TRUE
	)
)
