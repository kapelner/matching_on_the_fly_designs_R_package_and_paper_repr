#' Multivariate Log-Binomial Inference for Binary Responses
#'
#' @description
#' Fits a log-binomial model for binary (incidence) responses using treatment and
#' all recorded covariates as predictors. The treatment effect is reported on the
#' log-risk-ratio scale.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
DesignInferenceIncidMultiLogBinomial = R6::R6Class("DesignInferenceIncidMultiLogBinomial",
	inherit = DesignInferenceIncidUnivLogBinomial,
	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
