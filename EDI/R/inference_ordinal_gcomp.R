#' G-Computation Mean-Difference Inference for Ordinal Responses
#'
#' Fits a proportional-odds working model for an ordinal outcome using treatment
#' and, optionally, all recorded covariates, then estimates the marginal mean
#' difference by standardizing predicted mean ranks under all-treated and
#' all-control assignments over the empirical covariate distribution.
#'
#' @export
InferenceOrdinalGCompMeanDiff = R6::R6Class("InferenceOrdinalGCompMeanDiff",
	lock_objects = FALSE,
	inherit = InferenceOrdinalGCompAbstract,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = matrix(private$w, ncol = 1L)
				colnames(X) = "treatment"
			} else {
				X = cbind(treatment = private$w, X_cov)
			}
			X
		}
	)
)
