#' Inference based on Maximum Likelihood for KK designs
#'
#' @description
#' Inference for mean difference
#'
#'
#' @export
SeqDesignInferenceBaiAdjustedTKK21 = R6::R6Class("SeqDesignInferenceBaiAdjustedTKK21",
	inherit = SeqDesignInferenceBaiAdjustedT,
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
	#' @param convex_flag       A flag indicating whether the estimator should use a convex
	#'   combination of the Bai et al
	#'                                         matched pairs estimate with the reservoir estimate,
	#' or just the Bai et al estimate by its self.
	#'
	#' @examples
	#' set.seed(1)
	#' x_dat <- data.frame(
	#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
	#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
	#' )
	#' seq_des <- SeqDesignKK21$new(n = nrow(x_dat), response_type = "continuous", verbose = FALSE)
	#' for (i in seq_len(nrow(x_dat))) {
	#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
	#' }
	#' seq_des$add_all_subject_responses(c(1.2, 0.9, 1.5, 1.8, 2.1, 1.7, 2.6, 2.2))
	#' infer <- SeqDesignInferenceBaiAdjustedTKK21$new(seq_des, verbose = FALSE)
	#' infer
	#'
	initialize = function(seq_des_obj, num_cores = 1, verbose = TRUE, convex_flag = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose, convex_flag = convex_flag)
			assertNoCensoring(private$any_censoring)
		assert(checkClass(seq_des_obj, "SeqDesignKK21"), checkClass(seq_des_obj, "SeqDesignKK21stepwise"))
	}
	),

	private = list(
	distance = function(avg1, avg2){
		sum(private$seq_des_obj_priv_int$covariate_weights * (avg1 - avg2)^2)
	}
	)
)
