#' Cumulative Probit Inference for Ordinal Responses
#'
#' Ordinal probit (cumulative probit) model inference for ordinal responses.
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
#' infer <- InferenceOrdinalUniCumulProbitRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceOrdinalUniCumulProbitRegr = R6::R6Class("InferenceOrdinalUniCumulProbitRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj         A DesignSeqOneByOne object whose entire n subjects are assigned
		#'   and response y is recorded within.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{TRUE}.
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
		}
		),

	private = list(
		cumulative_probit_design_matrix = function(){
			Xmm = matrix(private$w, ncol = 1)
			colnames(Xmm) = c("treatment")
			Xmm
		},

		generate_mod = function(estimate_only = FALSE){
			if (estimate_only) {
				res = fast_ordinal_probit_regression_cpp(
					X = private$cumulative_probit_design_matrix(),
					y = as.numeric(private$y)
				)
				return(list(
					b = c(NA, res$b[1]),
					ssq_b_2 = NA_real_
				))
			}

			res = fast_ordinal_probit_regression_with_var_cpp(
				X = private$cumulative_probit_design_matrix(),
				y = as.numeric(private$y)
			)

			list(
				b = c(NA, res$b[1]),
				ssq_b_2 = res$ssq_b_2
			)
		}
	)
) # End of R6::R6Class
