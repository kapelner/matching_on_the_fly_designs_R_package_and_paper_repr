#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
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
#'   response_type = "proportion",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0.10, 0.25, 0.20, 0.40, 0.35, 0.55, 0.60, 0.75))
#' infer <- InferencePropMultiBetaRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferencePropMultiBetaRegr = R6::R6Class("InferencePropMultiBetaRegr",
	lock_objects = FALSE,
	inherit = InferencePropUniBetaRegr,
	public = list(

	),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			Xmm = private$create_design_matrix()
			# create_design_matrix is [Intercept, Treatment, Covariates]
			colnames(Xmm) = c("(Intercept)", "treatment", if(ncol(Xmm) > 2) paste0("x", 1:(ncol(Xmm)-2)) else NULL)
			y_beta = private$sanitize_beta_response(private$y)

			if (estimate_only) {
				res = fast_beta_regression(Xmm = Xmm, y = y_beta)
				# Ensure names are set for shared()
				names(res$b) = colnames(Xmm)
				return(list(
					b = res$b,
					ssq_b_2 = NA_real_
				))
			} else {
				res = fast_beta_regression_with_var(Xmm = Xmm, y = y_beta)
				# Ensure names are set for shared()
				names(res$b) = colnames(Xmm)
				return(list(
					b = res$b,
					ssq_b_2 = res$ssq_b_2
				))
			}
		}
	)
	)
