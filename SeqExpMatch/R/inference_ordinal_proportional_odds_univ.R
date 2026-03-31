#' Proportional Odds Inference based on Maximum Likelihood
#'
#' Proportional odds model (cumulative logit) inference for ordinal responses.
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
#' infer <- InferenceOrdinalUniPropOddsRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceOrdinalUniPropOddsRegr = R6::R6Class("InferenceOrdinalUniPropOddsRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects
		#'   are assigned and response y is recorded within.
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
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
		}
		),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			# Proportional odds model logit(P(Y <= k)) = alpha_k - beta * w
			# We use fast_ordinal_regression_with_var_cpp
			# Xmm should NOT include intercept as fast_ordinal_regression_cpp handles intercepts (alphas)
			Xmm = matrix(private$w, ncol = 1)
			colnames(Xmm) = c("treatment")
			
			if (estimate_only) {
				res = fast_ordinal_regression_cpp(X = Xmm, y = as.numeric(private$y))
				return(list(
					b = c(NA, res$b[1]), # Match the [2] indexing in shared()
					ssq_b_2 = NA_real_
				))
			}

			res = fast_ordinal_regression_with_var_cpp(X = Xmm, y = as.numeric(private$y))

			# Return in expected format
			# b[1] is the coefficient for treatment
			list(
				b = c(NA, res$b[1]), # Match the [2] indexing in shared()
				ssq_b_2 = res$ssq_b_2
			)
		}
	)
) # End of R6::R6Class
