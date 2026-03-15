#' Multivariate Hurdle Negative Binomial Regression Inference for Count Responses
#'
#' @description
#' Fits a hurdle negative binomial regression for count responses using the
#' treatment indicator and all recorded covariates in both the hurdle and count
#' components. The hurdle indicator model is fit on all subjects, and the
#' zero-truncated negative binomial count component is optimized in C++ on the
#' positive-count subjects. The reported treatment effect is the treatment
#' coefficient from the conditional count component, on the log-rate scale.
#'
#' @export
SeqDesignInferenceCountMultiHurdleNegBinRegr = R6::R6Class("SeqDesignInferenceCountMultiHurdleNegBinRegr",
	inherit = SeqDesignInferenceCountUnivHurdleNegBinRegr,
	private = list(
		predictors_df = function(){
			full_X = private$create_design_matrix()
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		}
	)
)
