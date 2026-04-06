#' Multivariate Linear Mixed Model Inference for KK Designs with Continuous Response
#'
#' Fits a linear mixed model using the \pkg{glmmTMB} fitter for continuous responses
#' under a KK matching-on-the-fly design using the treatment indicator and all
#' recorded covariates as fixed-effect predictors. A Gaussian identity-link working
#' model is used. The matched-pair strata enter the model as a subject-level random
#' intercept \code{(1 | group_id)}, which accounts for within-pair correlation.
#' Reservoir subjects each receive a unique singleton group. The treatment estimate
#' is the conditional mean difference; inference uses the Wald Z-statistic from the
#' mixed-model fixed-effects table.
#'
#' @details
#' This class requires the \pkg{glmmTMB} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{glmmTMB} before using this class.
#'
#' @export
InferenceContinMultiKKGLMM = R6::R6Class("InferenceContinMultiKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(




	),
	private = list(
		glmm_response_type = function() "continuous",
		glmm_family        = function() gaussian(link = "identity")
	)
)
