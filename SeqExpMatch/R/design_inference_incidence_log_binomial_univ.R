#' Univariate Log-Binomial Inference for Binary Responses
#'
#' @description
#' Fits a log-binomial model for binary (incidence) responses using only the
#' treatment indicator. The treatment effect is reported on the log-risk-ratio
#' scale.
#'
#' @details
#' The model uses a Bernoulli likelihood with log link, so the treatment
#' coefficient targets a log risk ratio directly. Fitting uses a constrained
#' IRLS routine implemented in C++ to respect the log-binomial mean constraint
#' \eqn{\mu_i < 1}.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
DesignInferenceIncidUnivLogBinomial = R6::R6Class("DesignInferenceIncidUnivLogBinomial",
	inherit = DesignInferenceIncidLogBinomialAbstract,
	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		}
	)
)
