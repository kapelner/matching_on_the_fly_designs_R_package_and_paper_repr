#' Univariate GLMM Inference for KK Designs with Binary Response
#'
#' Fits a Generalized Linear Mixed Model (GLMM) using the \pkg{glmmTMB} fitter for binary
#' (incidence) responses under a KK matching-on-the-fly design using only the
#' treatment indicator as a fixed-effect predictor (intercept + treatment). The
#' matched-pair strata enter the model as a subject-level random intercept
#' \code{(1 | group_id)}, which accounts for within-pair correlation. Reservoir
#' subjects each receive a unique singleton group. The treatment estimate is the
#' conditional (subject-specific) log-odds ratio; inference uses the Wald Z-statistic
#' from the GLMM fixed-effects table.
#'
#' @details
#' This class requires the \pkg{glmmTMB} package, which is listed in Suggests
#' and is not installed automatically with \pkg{SeqExpMatch}.
#' Install \pkg{glmmTMB} before using this class.
#'
#' @export
InferenceIncidUnivKKGLMM = R6::R6Class("InferenceIncidUnivKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(




	),
	private = list(
		glmm_response_type  = function() "incidence",
		glmm_family         = function() binomial(link = "logit"),
		glmm_predictors_df  = function() data.frame(w = private$w)
	)
)
