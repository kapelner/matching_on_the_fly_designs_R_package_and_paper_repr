#' GLMM Inference for KK Designs with Binary Response
#'
#' Fits a Generalized Linear Mixed Model (GLMM) using the \pkg{glmmTMB} fitter for
#' binary (incidence) responses under a KK matching-on-the-fly design using the
#' treatment indicator and, optionally, all recorded covariates as fixed-effect
#' predictors.
#'
#' @details
#' This class requires the \pkg{glmmTMB} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{glmmTMB} before using this class.
#'
#' @export
InferenceIncidKKGLMM = R6::R6Class("InferenceIncidKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(
	),
	private = list(
		glmm_response_type  = function() "incidence",
		glmm_family         = function() stats::binomial(link = "logit")
	)
)
