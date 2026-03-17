#' Binomial Identity-Link Inference for Binary Responses
#'
#' @description
#' Internal base class for binomial identity-link regression inference on
#' incidence outcomes. The treatment effect is reported on the risk-difference
#' scale.
#'
#' @keywords internal
#' @noRd
SeqDesignInferenceIncidBinomialIdentityAbstract = R6::R6Class("SeqDesignInferenceIncidBinomialIdentityAbstract",
	inherit = SeqDesignInferenceIncidConstrainedBinomialAbstract,
	private = list(
		fit_constrained_binomial = function(X_fit, j_treat){
			fast_identity_binomial_regression_with_var_cpp(X_fit, as.numeric(private$y), j = j_treat)
		},
		method_label = function() "Binomial identity-link regression"
	)
)
