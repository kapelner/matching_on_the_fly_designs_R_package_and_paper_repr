#' Multivariate KK Hurdle Poisson Combined-Likelihood Inference for Count Responses
#'
#' @description
#' Fits a KK hurdle-Poisson combined-likelihood model for count responses using
#' the treatment indicator and all recorded covariates in both the hurdle and
#' positive-count components. Matched pairs contribute pair-specific random
#' intercepts, while reservoir subjects contribute through a common fixed
#' reservoir intercept. The reported treatment effect is the treatment
#' coefficient from the positive-count component on the log-rate scale.
#'
#' @export
InferenceCountMultiKKHurdlePoissonCombinedLikelihood = R6::R6Class("InferenceCountMultiKKHurdlePoissonCombinedLikelihood",
	inherit = InferenceAbstractKKHurdlePoissonCombinedLikelihood,
	private = list(
		include_covariates = function() TRUE,

		predictors_df = function(){
			full_X = private$create_design_matrix()
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		}
	)
)
