#' G-Computation Risk-Difference Inference for Binary Responses
#'
#' Fits a logistic working model for an incidence outcome using treatment and,
#' optionally, all recorded covariates, then estimates the marginal risk
#' difference by standardizing predicted risks under all-treated and all-control
#' assignments over the empirical covariate distribution.
#'
#' @export
InferenceIncidGCompRiskDiff = R6::R6Class("InferenceIncidGCompRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceIncidGCompAbstract,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			if (private$include_covariates) {
				private$create_design_matrix()
			} else {
				cbind("(Intercept)" = 1, treatment = private$w)
			}
		},

		get_estimand_type = function() "RD"
	)
)

#' G-Computation Risk-Ratio Inference for Binary Responses
#'
#' Fits a logistic working model for an incidence outcome using treatment and,
#' optionally, all recorded covariates, then estimates the marginal risk ratio
#' by standardizing predicted risks under all-treated and all-control assignments
#' over the empirical covariate distribution.
#'
#' @export
InferenceIncidGCompRiskRatio = R6::R6Class("InferenceIncidGCompRiskRatio",
	lock_objects = FALSE,
	inherit = InferenceIncidGCompAbstract,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			if (private$include_covariates) {
				private$create_design_matrix()
			} else {
				cbind("(Intercept)" = 1, treatment = private$w)
			}
		},

		get_estimand_type = function() "RR"
	)
)
