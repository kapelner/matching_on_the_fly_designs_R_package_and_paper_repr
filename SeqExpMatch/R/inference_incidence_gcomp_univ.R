#' Univariate G-Computation Risk-Difference Inference for Binary Responses
#'
#' Fits a univariate logistic working model for an incidence outcome and then
#' estimates the marginal risk difference by standardizing predicted risks under
#' all-treated and all-control assignments over the empirical covariate
#' distribution.
#'
#' @export
InferenceIncidUnivGCompRiskDiff = R6::R6Class("InferenceIncidUnivGCompRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceIncidGCompAbstract,
	public = list(





	),

	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		},

		get_estimand_type = function() "RD"
	)
)
InferenceIncidUnivGCompRiskRatio = R6::R6Class("InferenceIncidUnivGCompRiskRatio",
	lock_objects = FALSE,
	inherit = InferenceIncidGCompAbstract,
	public = list(





	),

	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		},

		get_estimand_type = function() "RR"
	)
)
