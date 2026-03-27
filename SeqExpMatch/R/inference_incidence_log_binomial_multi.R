#' Multivariate Log-Binomial Inference for Binary Responses
#'
#' @description
#' Fits a log-binomial model for binary (incidence) responses using treatment and
#' all recorded covariates as predictors. The treatment effect is reported on the
#' log-risk-ratio scale.
#'
#' @export
InferenceIncidMultiLogBinomial = R6::R6Class("InferenceIncidMultiLogBinomial",
	inherit = InferenceIncidUnivLogBinomial,
	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
