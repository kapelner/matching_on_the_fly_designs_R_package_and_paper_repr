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
		#' @param seq_des_obj             A SeqDesign object whose entire n subjects are assigned and
		#'   whose continuous response y is recorded.
		#' @param method                  Robust-regression fitting method for \code{MASS::rlm}; one
		#'   of \code{"M"} or \code{"MM"}. The default is \code{"MM"}.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#' 							and bootstrap resampling.
		#' @param verbose                 A flag indicating whether messages should be displayed to
		#'   the user. Default is \code{FALSE}.
		#' @examples
		#' set.seed(1)
		#' x_dat <- data.frame(
		#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
		#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
		#' )
		#' seq_des <- SeqDesignCRD$new(n = nrow(x_dat), response_type = "continuous", verbose = FALSE)
		#' for (i in seq_len(nrow(x_dat))) {
		#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
		#' }
		#' seq_des$add_all_subject_responses(c(1.2, 0.9, 1.5, 1.8, 2.1, 1.7, 2.6, 2.2))
		#' infer <- SeqDesignInferenceContinMultiRobustRegr$new(seq_des, verbose = FALSE)
		#' infer
		#'
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
