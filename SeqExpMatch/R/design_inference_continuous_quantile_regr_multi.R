#' Multivariate Quantile Regression Inference for Continuous Responses
#'
#' @description
#' Fits a quantile regression for continuous responses using the treatment
#' indicator and all recorded covariates as predictors. The treatment effect is
#' reported on the response scale at quantile \code{tau}; by default
#' \code{tau = 0.5}, so this is median regression.
#'
#' @export
DesignInferenceContinMultiQuantileRegr = R6::R6Class("DesignInferenceContinMultiQuantileRegr",
	inherit = DesignInferenceContinUnivQuantileRegr,
	public = list(

		#' @description
		#' Initialize a multivariate quantile-regression inference object for a completed
		#' design with a continuous response.
		#' @param des_obj A completed \code{SeqDesign} object with a continuous response.
		#' @param tau The quantile level for regression, strictly between 0 and 1. The default is
		#'   \code{tau = 0.5}.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		#' @examples
		#' set.seed(1)
		#' x_dat <- data.frame(
		#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
		#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
		#' )
		#' seq_des <- SeqDesignBernoulli$new(n = nrow(x_dat), response_type = "continuous", verbose = FALSE)
		#' for (i in seq_len(nrow(x_dat))) {
		#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
		#' }
		#' seq_des$add_all_subject_responses(c(1.2, 0.9, 1.5, 1.8, 2.1, 1.7, 2.6, 2.2))
		#' infer <- DesignInferenceContinMultiQuantileRegr$new(seq_des, verbose = FALSE)
		#' infer
		#'
		initialize = function(des_obj, tau = 0.5, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, tau = tau, num_cores = num_cores, verbose = verbose)
		}
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
