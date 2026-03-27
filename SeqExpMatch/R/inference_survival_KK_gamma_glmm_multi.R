#' Multivariate GLMM Inference for KK Designs with Survival Response (no censoring)
#'
#' @description
#' Fits a Generalized Linear Mixed Model (GLMM) using the \pkg{glmmTMB} fitter for uncensored
#' survival (time-to-event) responses under a KK matching-on-the-fly design using the
#' treatment indicator and all recorded covariates as fixed-effect predictors. A
#' Gamma log-link working model is used, which is appropriate for positive continuous
#' survival times. The matched-pair strata enter the model as a subject-level random
#' intercept \code{(1 | group_id)}, which accounts for within-pair correlation.
#' Reservoir subjects each receive a unique singleton group. The treatment estimate is
#' the conditional log ratio of mean survival times (log-MTR); inference uses the Wald
#' Z-statistic from the GLMM fixed-effects table.
#'
#' @details
#' Censored observations are not supported; an error is raised at initialization if
#' any censoring is detected. This class requires the \pkg{glmmTMB} package, which is
#' listed under \code{Suggests} and is not installed automatically with
#' \pkg{SeqExpMatch}. Install \pkg{glmmTMB} manually
#' before using this class.
#'
#' @export
InferenceSurvivalMultiKKGammaGLMM = R6::R6Class("InferenceSurvivalMultiKKGammaGLMM",
	inherit = InferenceAbstractKKGLMM,
	public = list(




	),
	private = list(
		glmm_response_type = function() "survival",
		glmm_family        = function() Gamma(link = "log")
	)
)
