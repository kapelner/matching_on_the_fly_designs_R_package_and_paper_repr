#' G-Computation Mean-Difference Inference for Proportion Responses
#'
#' Fits a fractional-logit working model for a proportion outcome using treatment
#' and, optionally, all recorded covariates, then estimates the marginal mean
#' difference by standardizing predicted mean proportions under all-treated and
#' all-control assignments over the empirical covariate distribution.
#'
#' @export
InferencePropGCompMeanDiff = R6::R6Class("InferencePropGCompMeanDiff",
	lock_objects = FALSE,
	inherit = InferencePropGCompAbstract,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			if (private$include_covariates) {
				private$create_design_matrix()
			} else {
				cbind("(Intercept)" = 1, treatment = private$w)
			}
		}
	)
)
