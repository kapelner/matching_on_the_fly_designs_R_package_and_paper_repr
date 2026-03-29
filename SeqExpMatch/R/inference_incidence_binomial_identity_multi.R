#' Multivariate Binomial Identity-Link Regression for Binary Responses
#'
#' Fits a binomial GLM with identity link for binary (incidence) responses using
#' treatment and all recorded covariates as predictors. The treatment effect is
#' reported on the risk-difference scale.
#'
#' @export
InferenceIncidMultiBinomialIdentityRiskDiff = R6::R6Class("InferenceIncidMultiBinomialIdentityRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceIncidUnivBinomialIdentityRiskDiff,
	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
