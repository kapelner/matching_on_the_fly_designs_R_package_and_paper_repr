#' Multivariate Robust Regression Inference for Continuous Responses
#'
#' @description
#' Fits a robust linear regression via \code{MASS::rlm} for continuous responses
#' using the treatment indicator and all recorded covariates as predictors.
#' This provides a Huber/MM-style robustness upgrade over ordinary least squares
#' when outcomes are heavy-tailed or outlier-prone. Inference is based on the
#' coefficient table returned by \code{summary.rlm()}.
#'
#' @export
SeqDesignInferenceContinMultiRobustRegr = R6::R6Class("SeqDesignInferenceContinMultiRobustRegr",
	inherit = SeqDesignInferenceContinUnivRobustRegr,
	public = list(

		#' @description
		#' Initialize a multivariate robust-regression inference object for a completed design
		#' with a continuous response.
		#' @param	seq_des_obj		A SeqDesign object whose entire n subjects are assigned and whose continuous response y is recorded.
		#' @param	method			Robust-regression fitting method for \code{MASS::rlm}; one of \code{"M"} or \code{"MM"}. The default is \code{"MM"}.
		#' @param	num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling.
		#' @param	verbose			A flag indicating whether messages should be displayed to the user. Default is \code{FALSE}.
		initialize = function(seq_des_obj, method = "MM", num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, method = method, num_cores = num_cores, verbose = verbose)
		}
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
