#' Sequential Design Inference for Count Response Types using Multivariate Negative Binomial
#' Regression
#'
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
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
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
	lock_objects = FALSE,
	inherit = InferenceCountUnivNegBinRegr,
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
						res = fast_negbin_regression(Xmm = X_fit, y = private$y)
						list(b = res$b, ssq_b_2 = NA_real_, j_treat = j_treat)
					} else {
						res = fast_negbin_regression_with_var(Xmm = X_fit, y = private$y)
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
