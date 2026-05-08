#' G-Computation Mean-Difference Inference for Proportion Responses
#'
#' Fits a fractional-logit working model for a proportion outcome using treatment
#' and, optionally, all recorded covariates, then estimates the marginal mean
#' difference by standardizing predicted mean proportions under all-treated and
#' all-control assignments over the empirical covariate distribution.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'proportion')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferencePropGCompMeanDiff$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferencePropGCompMeanDiff = R6::R6Class("InferencePropGCompMeanDiff",
	lock_objects = FALSE,
	inherit = InferencePropGCompAbstract,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
