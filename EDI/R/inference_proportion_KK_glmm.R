#' GLMM Inference for KK Designs with Proportion Response
#'
#' Fits a Generalized Linear Mixed Model (GLMM) using the \pkg{glmmTMB} fitter for
#' proportion (continuous values in (0, 1)) responses under a KK
#' matching-on-the-fly design using the treatment indicator and, optionally, all
#' recorded covariates as fixed-effect predictors.
#'
#' @details
#' This class requires the \pkg{glmmTMB} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{glmmTMB} before using this class.
#'
#' @export
InferencePropKKGLMM = R6::R6Class("InferencePropKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(
	),
	private = list(
		glmm_response_type  = function() "proportion",
		glmm_family         = function() stats::binomial(link = "logit")
	)
)
