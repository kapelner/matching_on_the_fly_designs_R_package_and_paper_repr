#' Continuation-Ratio Logit Inference based on Maximum Likelihood
#'
#' Continuation-ratio logit model inference for ordinal responses. This model
#' conditions on "having reached category k" before comparing category k vs >k.
#' It is particularly useful when categories represent a progression.
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
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalContRatioRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceOrdinalContRatioRegr = R6::R6Class("InferenceOrdinalContRatioRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj         A DesignSeqOneByOne object whose entire n subjects are assigned
		#'   and response y is recorded within.
		#' @param verbose                 A flag indicating whether messages should be displayed
		#'   to the user. Default is \code{TRUE}.
		#' @param harden Whether to apply robustness measures.
		initialize = function(des_obj,  verbose = FALSE, harden = TRUE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, verbose, harden)
			assertNoCensoring(private$any_censoring)
		}
		),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			# Continuation-ratio model logit(P(Y = k | Y >= k)) = alpha_k + beta * w
			# We use fast_continuation_ratio_regression_with_var_cpp
			Xmm = matrix(private$w, ncol = 1)
			full_names = c("treatment")
			colnames(Xmm) = full_names[seq_len(ncol(Xmm))]
			
			if (estimate_only) {
				res = fast_continuation_ratio_regression_cpp(X = Xmm, y = as.numeric(private$y))
				return(list(
					b = c(NA, res$b[1]), # Match the [2] indexing in shared()
					ssq_b_2 = NA_real_
				))
			}

			res = fast_continuation_ratio_regression_with_var_cpp(X = Xmm, y = as.numeric(private$y))

			# Return in expected format
			list(
				b = c(NA, res$b[1]), # Match the [2] indexing in shared()
				ssq_b_2 = res$ssq_b_2
			)
		}
	)
) # End of R6::R6Class
