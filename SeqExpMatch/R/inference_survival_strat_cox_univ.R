#' Univariate Stratified Cox PH Regression Inference for Survival Responses
#'
#' @description
#' Fits an all-subject stratified Cox proportional hazards model for survival
#' responses using only the treatment indicator in the linear predictor. The strata
#' are chosen automatically from low-cardinality observed covariates. If no suitable
#' stratification covariates are available, the method falls back to the standard
#' univariate Cox PH model.
#'
#' @export
SeqDesignInferenceSurvivalUniStratCoxPHRegr = R6::R6Class("SeqDesignInferenceSurvivalUniStratCoxPHRegr",
	inherit = SeqDesignInferenceSurvivalStratCoxPHAbstract,
	private = list(
		include_covariates = function() FALSE
	)
)
