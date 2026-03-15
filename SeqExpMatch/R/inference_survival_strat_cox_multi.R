#' Multivariate Stratified Cox PH Regression Inference for Survival Responses
#'
#' @description
#' Fits an all-subject stratified Cox proportional hazards model for survival
#' responses using the treatment indicator and the recorded covariates in the
#' linear predictor. Low-cardinality covariates are reserved for automatic
#' stratification and the remaining covariates are included as regression terms.
#' If no suitable stratification covariates are available, the method falls back to
#' the standard multivariate Cox PH model.
#'
#' @export
SeqDesignInferenceSurvivalMultiStratCoxPHRegr = R6::R6Class("SeqDesignInferenceSurvivalMultiStratCoxPHRegr",
	inherit = SeqDesignInferenceSurvivalStratCoxPHAbstract,
	private = list(
		include_covariates = function() TRUE
	)
)
