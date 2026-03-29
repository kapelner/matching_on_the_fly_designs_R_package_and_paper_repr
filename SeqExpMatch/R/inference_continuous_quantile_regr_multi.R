#' Multivariate Quantile Regression Inference for Continuous Responses
#'
#' Fits a quantile regression for continuous responses using the treatment
#' indicator and all recorded covariates as predictors. The treatment effect is
#' reported on the response scale at quantile \code{tau}; by default
#' \code{tau = 0.5}, so this is median regression.
#'
#' @export
InferenceContinMultiQuantileRegr = R6::R6Class("InferenceContinMultiQuantileRegr",
	lock_objects = FALSE,
	inherit = InferenceContinUnivQuantileRegr,
	public = list(

	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
