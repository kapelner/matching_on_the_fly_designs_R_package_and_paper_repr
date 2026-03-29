#' Cumulative Complementary Log-Log Inference based on Maximum Likelihood
#'
#' @description
#' Cumulative complementary log-log model inference for ordinal responses.
#' This model is also known as the discrete-time proportional hazards model.
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
#'   response_type = "ordinal",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalUniCLLRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceOrdinalUniCLLRegr = R6::R6Class("InferenceOrdinalUniCLLRegr",
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
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
		}
		),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			# Cumulative CLL model cloglog(P(Y <= k)) = alpha_k + beta * w
			# We use fast_ordinal_cloglog_regression_with_var_cpp
			Xmm = matrix(private$w, ncol = 1)
			colnames(Xmm) = c("treatment")
			
			if (estimate_only) {
				res = fast_ordinal_cloglog_regression_cpp(X = Xmm, y = as.numeric(private$y))
				return(list(
					b = c(NA, res$b[1]), # Match the [2] indexing in shared()
					ssq_b_2 = NA_real_
				))
			}

			res = fast_ordinal_cloglog_regression_with_var_cpp(X = Xmm, y = as.numeric(private$y))

			# Return in expected format
			list(
				b = c(NA, res$b[1]), # Match the [2] indexing in shared()
				ssq_b_2 = res$ssq_b_2
			)
		}
	)
) # End of R6::R6Class
