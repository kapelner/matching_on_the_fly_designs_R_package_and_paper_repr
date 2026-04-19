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
			# Fast path for bootstrap/jackknife: reuse the column structure from the
			# original fit instead of re-running the expensive create_design_matrix()
			# (drop_highly_correlated_cols + drop_linearly_dependent_cols) each time.
			if (estimate_only && !is.null(private$best_Xmm_colnames)) {
				X_data = private$get_X()
				Xmm_cols = private$best_Xmm_colnames
				if (length(Xmm_cols) == 0L) {
					Xmm = cbind(1, private$w)
				} else {
					X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
					Xmm = cbind(1, treatment = private$w, X_cov)
				}
				res = fast_logistic_regression_cpp(X = Xmm, y = as.numeric(private$y))
				if (is.null(res) || !is.finite(res$b[2])) return(list(b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_))
				return(list(b = res$b, ssq_b_2 = NA_real_))
			}

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
