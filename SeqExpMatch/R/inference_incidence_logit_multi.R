#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#'
#'
#' @export
InferenceIncidMultiLogRegr = R6::R6Class("InferenceIncidMultiLogRegr",
	inherit = InferenceIncidUnivLogRegr,
	public = list(

	),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			if (estimate_only) {
				res = fast_logistic_regression(private$create_design_matrix(), private$y)
				list(b = res$b, ssq_b_2 = NA_real_)
			} else {
				fast_logistic_regression_with_var(private$create_design_matrix(), private$y)
			}
		}
	)
)
