#' Multivariate Robust Regression Inference for Continuous Responses
#'
#' Fits a robust linear regression via \code{MASS::rlm} for continuous responses
#' using the treatment indicator and all recorded covariates as predictors.
#' This provides a Huber/MM-style robustness upgrade over ordinary least squares
#' when outcomes are heavy-tailed or outlier-prone. Inference is based on the
#' coefficient table returned by \code{summary.rlm()}.
#'
#' @export
InferenceContinMultiRobustRegr = R6::R6Class("InferenceContinMultiRobustRegr",
	lock_objects = FALSE,
	inherit = InferenceContinUnivRobustRegr,
	public = list(

	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
