#' Log-Binomial Inference for Binary Responses
#'
#' @description
#' Internal base class for log-binomial regression inference on incidence
#' outcomes. The treatment effect is reported on the log-risk-ratio scale.
#'
#' @keywords internal
#' @noRd
DesignInferenceIncidLogBinomialAbstract = R6::R6Class("DesignInferenceIncidLogBinomialAbstract",
	inherit = DesignInferenceAsymp,
	private = list(
		fit_constrained_binomial = function(X_fit, j_treat){
			fast_log_binomial_regression_with_var_cpp(X_fit, as.numeric(private$y), j = j_treat)
		},
		method_label = function() "Log-binomial estimator"
	)
)
