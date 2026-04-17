#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#'
#'
#' @export
InferenceIncidMultiLogRegr = R6::R6Class("InferenceIncidMultiLogRegr",
	lock_objects = FALSE,
	inherit = InferenceIncidUnivLogRegr,
	public = list(

	),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			X_full = private$create_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit, j_treat){
					if (estimate_only) {
						res = fast_logistic_regression_cpp(X_fit, private$y)
						list(b = res$b, ssq_b_2 = NA_real_, j_treat = j_treat)
					} else {
						res = fast_logistic_regression_with_var_cpp(X_fit, private$y, j = j_treat)
						res$j_treat = j_treat
						res
					}
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = mod$j_treat
					if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)

			if (!is.null(attempt$fit)){
				private$best_Xmm_colnames = setdiff(colnames(attempt$X_fit), c("(Intercept)", "treatment"))
			}
			attempt$fit
		}
	)
)
