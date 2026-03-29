#' Multivariate G-Computation Mean-Difference Inference for Proportion Responses
#'
#' Fits a multivariate fractional-logit working model for a proportion outcome
#' and then estimates the marginal mean difference by standardizing predicted
#' mean proportions under all-treated and all-control assignments over the
#' empirical covariate distribution.
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
#' infer <- InferencePropMultiGCompMeanDiff$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferencePropMultiGCompMeanDiff = R6::R6Class("InferencePropMultiGCompMeanDiff",
	lock_objects = FALSE,
	inherit = InferencePropUniGCompMeanDiff,
	public = list(

	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
