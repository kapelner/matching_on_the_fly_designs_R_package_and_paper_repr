#' Sequential Design Inference for Count Response Types using Multivariate Negative Binomial
#' Regression
#'
#' @description
#' The methods that support confidence intervals and testing for
#' count response Types using Multivariate Negative Binomial Regression
#'
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignCRD$new(n = nrow(x_dat), response_type = "count", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(c(0, 1, 1, 2, 2, 3, 3, 4))
#' infer <- SeqDesignInferenceCountMultiNegBinRegr$new(seq_des, verbose = FALSE)
#' infer
#'
SeqDesignInferenceCountMultiNegBinRegr = R6::R6Class("SeqDesignInferenceCountMultiNegBinRegr",
	inherit = SeqDesignInferenceCountUnivNegBinRegr,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param seq_des_obj A SeqDesign object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#'   The default is 1 for serial computation. For simple
		#'   estimators (e.g. mean difference and KK compound),
		#'   parallelization is achieved with zero-overhead C++ OpenMP.
		#'   For complex models (e.g. GLMs),
		#'   parallelization falls back to R's
		#'   \code{parallel::mclapply}, which incurs
		#'   session-forking overhead.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{TRUE}.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),

	private = list(
		generate_mod = function(){
			fast_negbin_regression_with_var(
				Xmm = private$create_design_matrix(),
				y = private$y
			)
		}
	)
)
