#' Multivariate Hurdle Poisson Regression Inference for Count Responses
#'
#' @description
#' Fits a hurdle Poisson regression for count responses using the treatment
#' indicator and all recorded covariates in both the count and hurdle components.
#' The reported treatment effect is the treatment coefficient from the conditional
#' truncated-Poisson count component, on the log-rate scale.
#'
#' @export
SeqDesignInferenceCountMultiHurdlePoissonRegr = R6::R6Class("SeqDesignInferenceCountMultiHurdlePoissonRegr",
	inherit = SeqDesignInferenceCountUnivHurdlePoissonRegr,
	private = list(
		predictors_df = function(){
			full_X = private$create_design_matrix()
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		}
	)
)
