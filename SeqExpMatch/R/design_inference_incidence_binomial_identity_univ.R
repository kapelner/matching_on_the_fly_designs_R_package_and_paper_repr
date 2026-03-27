#' Univariate Binomial Identity-Link Regression for Binary Responses
#'
#' @description
#' Fits a binomial GLM with identity link for binary (incidence) responses using
#' only the treatment indicator. The treatment effect is reported on the
#' risk-difference scale.
#'
#' @details
#' This model targets the adjusted risk difference directly through a Bernoulli
#' likelihood with identity link. Fitting uses a constrained IRLS routine
#' implemented in C++ to enforce fitted probabilities in \eqn{(0,1)}.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
DesignInferenceIncidUnivBinomialIdentityRiskDiff = R6::R6Class("DesignInferenceIncidUnivBinomialIdentityRiskDiff",
	inherit = DesignInferenceIncidBinomialIdentityAbstract,
	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		}
	)
)
