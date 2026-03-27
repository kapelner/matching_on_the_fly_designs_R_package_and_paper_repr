#' Univariate G-Computation Mean-Difference Inference for Proportion Responses
#'
#' @description
#' Fits a univariate fractional-logit working model for a proportion outcome and
#' then estimates the marginal mean difference by standardizing predicted mean
#' proportions under all-treated and all-control assignments over the empirical
#' covariate distribution.
#'
#' @export
InferencePropUniGCompMeanDiff = R6::R6Class("InferencePropUniGCompMeanDiff",
	inherit = InferencePropGCompAbstract,
	public = list(





	),

	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		}
	)
)
