#' Multivariate Zero-Inflated Poisson Regression Inference for Count Responses
#'
#' @description
#' Fits a zero-inflated Poisson regression for count responses using the treatment
#' indicator and all recorded covariates in both the count and zero-inflation
#' components. The reported treatment effect is the treatment coefficient from the
#' conditional Poisson count component, on the log-rate scale.
#'
#' @export
SeqDesignInferenceCountMultiZeroInflatedPoissonRegr = R6::R6Class("SeqDesignInferenceCountMultiZeroInflatedPoissonRegr",
	inherit = SeqDesignInferenceCountUnivZeroInflatedPoissonRegr,
	private = list(
		predictors_df = function(){
			full_X = private$create_design_matrix()
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		}
	)
)
