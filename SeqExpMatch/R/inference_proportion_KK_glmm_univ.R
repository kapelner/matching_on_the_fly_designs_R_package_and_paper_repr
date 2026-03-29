#' Univariate GLMM Inference for KK Designs with Proportion Response
#'
#' Fits a Generalized Linear Mixed Model (GLMM) using the \pkg{glmmTMB} fitter for proportion
#' (continuous values in (0, 1)) responses under a KK matching-on-the-fly design
#' using only the treatment indicator as a fixed-effect predictor (intercept +
#' treatment). A binomial logit working model is used; \code{glmmTMB} accepts continuous
#' proportions with this family (non-integer-success warnings are suppressed
#' internally). The matched-pair strata enter the model as a subject-level random
#' intercept \code{(1 | group_id)}, which accounts for within-pair correlation.
#' Reservoir subjects each receive a unique singleton group. The treatment estimate
#' is the conditional log-odds ratio of the mean proportion; inference uses the Wald
#' Z-statistic from the GLMM fixed-effects table.
#'
#' @details
#' This class requires the \pkg{glmmTMB} package, which is listed in Suggests
#' and is not installed automatically with \pkg{SeqExpMatch}.
#' Install \pkg{glmmTMB} before using this class.
#'
#' @export
InferencePropUnivKKGLMM = R6::R6Class("InferencePropUnivKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(




	),
	private = list(
		glmm_response_type  = function() "proportion",
		glmm_family         = function() binomial(link = "logit"),
		glmm_predictors_df  = function() data.frame(w = private$w)
	)
)
