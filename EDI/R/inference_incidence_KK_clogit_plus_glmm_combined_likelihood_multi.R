#' Multivariate Conditional Logistic Plus GLMM Combined-Likelihood Inference for KK Designs
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs. The treatment effect is estimated by the conditional-logistic
#' component; the GLMM component includes only the intercept and covariates.
#'
#' @details
#' The initializer accepts \code{max_abs_reasonable_coef} (default \code{1e4})
#' to bound finite treatment estimates and standard errors before declaring a fit
#' non-estimable, and \code{max_abs_log_sigma} (default \code{8}) to bound the
#' concordant-pair GLMM log random-intercept standard deviation.
#'
#' @export
InferenceIncidMultiKKClogitPlusGLMMCombinedLikelihood = R6::R6Class("InferenceIncidMultiKKClogitPlusGLMMCombinedLikelihood",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitPlusGLMM,
	public = list(),
	private = list(
		combine_reservoir_into_glmm = function() TRUE
	)
)
