#' Multivariate Quantile Regression Inference for Continuous Responses
#'
#' @description
#' Fits a quantile regression for continuous responses using the treatment
#' indicator and all recorded covariates as predictors. The treatment effect is
#' reported on the response scale at quantile \code{tau}; by default
#' \code{tau = 0.5}, so this is median regression.
#'
#' @export
SeqDesignInferenceContinMultiQuantileRegr = R6::R6Class("SeqDesignInferenceContinMultiQuantileRegr",
	inherit = SeqDesignInferenceContinUnivQuantileRegr,
	public = list(

		#' @description
		#' Initialize a multivariate quantile-regression inference object for a completed
		#' design with a continuous response.
		#' @param seq_des_obj A completed \code{SeqDesign} object with a continuous response.
		#' @param tau The quantile level for regression, strictly between 0 and 1. The default is \code{tau = 0.5}.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, tau = 0.5, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, tau = tau, num_cores = num_cores, verbose = verbose)
		}
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
