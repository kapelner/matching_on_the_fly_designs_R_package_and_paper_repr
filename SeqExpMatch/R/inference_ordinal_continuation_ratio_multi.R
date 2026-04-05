#' Multivariate Continuation-Ratio Logit Inference for Ordinal Responses
#'
#' Continuation-ratio logit model inference for ordinal responses with
#' baseline covariate adjustment.
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
#' infer <- InferenceOrdinalMultiContRatioRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalMultiContRatioRegr = R6::R6Class(
	"InferenceOrdinalMultiContRatioRegr",
	lock_objects = FALSE,
	inherit = InferenceOrdinalContRatioRegr,
	public = list(
	),

	private = list(
		contin_ratio_design_matrix = function(){
			X_full = cbind(private$w, private$get_X())
			qr_X = qr(X_full)
			if (qr_X$rank < ncol(X_full)){
				keep = qr_X$pivot[seq_len(qr_X$rank)]
				if (!(1L %in% keep)) keep[qr_X$rank] = 1L
				keep = sort(keep)
				X_full = X_full[, keep, drop = FALSE]
			}
			colnames(X_full)[1] = "treatment"
			X_full
		},

		generate_mod = function(estimate_only = FALSE){
			if (estimate_only) {
				res = fast_continuation_ratio_regression_cpp(
					X = private$contin_ratio_design_matrix(),
					y = as.numeric(private$y)
				)
				return(list(
					b = c(NA, res$b), # Match the [2] indexing in shared()
					ssq_b_2 = NA_real_
				))
			}

			res = fast_continuation_ratio_regression_with_var_cpp(
				X = private$contin_ratio_design_matrix(),
				y = as.numeric(private$y)
			)
			list(
				b = c(NA, res$b), # Match the [2] indexing in shared()
				ssq_b_2 = res$ssq_b_2
			)
		}
	)
)
