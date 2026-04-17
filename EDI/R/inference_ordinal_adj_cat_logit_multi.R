#' Multivariate Adjacent-Category Logit Inference for Ordinal Responses
#'
#' Adjacent-category logit inference for ordinal responses with treatment and
#' observed covariates entering linearly into the adjacent-category logits.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$new(n = nrow(x_dat), response_type = "ordinal",
#'   verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalMultiAdjCatLogitRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalMultiAdjCatLogitRegr = R6::R6Class(
	"InferenceOrdinalMultiAdjCatLogitRegr",
	lock_objects = FALSE,
	inherit = InferenceOrdinalUniAdjCatLogitRegr,
	public = list(

	),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			X_full = cbind(treatment = private$w, private$get_X())
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					if (estimate_only) {
						res = fast_adjacent_category_logit_cpp(
							X = X_fit,
							y = as.numeric(private$y)
						)
						return(list(
							b = c(NA, res$b),
							ssq_b_2 = NA_real_,
							converged = res$converged
						))
					}

					res = fast_adjacent_category_logit_with_var_cpp(
						X = X_fit,
						y = as.numeric(private$y)
					)
					list(
						b = c(NA, res$b),
						ssq_b_2 = res$ssq_b_1,
						converged = res$converged
					)
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
					if (!is.null(mod$converged) && !mod$converged) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$best_Xmm_colnames = setdiff(colnames(attempt$X_fit), "treatment")
			}
			attempt$fit
		}
	)
)
