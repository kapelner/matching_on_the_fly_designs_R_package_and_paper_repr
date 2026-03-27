#' Multivariate Binomial Identity-Link Regression for Binary Responses
#'
#' @description
#' Fits a binomial GLM with identity link for binary (incidence) responses using
#' treatment and all recorded covariates as predictors. The treatment effect is
#' reported on the risk-difference scale.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
DesignInferenceIncidMultiBinomialIdentityRiskDiff = R6::R6Class("DesignInferenceIncidMultiBinomialIdentityRiskDiff",
	inherit = DesignInferenceIncidUnivBinomialIdentityRiskDiff,
	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
