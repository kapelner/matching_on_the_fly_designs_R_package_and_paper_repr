#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#'
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "proportion",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0.10, 0.25, 0.20, 0.40, 0.35, 0.55, 0.60, 0.75))
#' infer <- InferencePropUniBetaRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferencePropUniBetaRegr = R6::R6Class("InferencePropUniBetaRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects
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
		#' @param make_fork_cluster Whether to use a fork cluster for parallelization.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "proportion")
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
			assertNumeric(private$y, any.missing = FALSE, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
		}
		),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			Xmm = cbind(1, private$w)
			colnames(Xmm) = c("(Intercept)", "treatment")
			
			if (estimate_only) {
				res = fast_beta_regression(Xmm = Xmm, y = private$y)
				# Ensure names are set for shared()
				names(res$b) = colnames(Xmm)
				return(list(
					b = res$b,
					ssq_b_2 = NA_real_
				))
			} else {
				res = fast_beta_regression_with_var(Xmm = Xmm, y = private$y)
				# Ensure names are set for shared()
				names(res$b) = colnames(Xmm)
				return(list(
					b = res$b,
					ssq_b_2 = res$ssq_b_2
				))
			}
		}
	)
) # End of R6::R6Class
