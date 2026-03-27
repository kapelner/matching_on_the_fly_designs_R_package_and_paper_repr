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
#' seq_des <- DesignSeqOneByOneBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "count",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 1, 2, 2, 3, 3, 4))
#' infer <- InferenceCountMultiNegBinRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceCountMultiNegBinRegr = R6::R6Class("InferenceCountMultiNegBinRegr",
	inherit = InferenceCountUnivNegBinRegr,
	public = list(

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
