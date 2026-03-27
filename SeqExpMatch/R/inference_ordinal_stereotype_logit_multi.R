#' Multivariate Stereotype Logistic Inference for Ordinal Responses
#'
#' @description
#' Stereotype logistic inference for ordinal responses with treatment and observed
#' covariates entering the linear predictor.
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
#'   response_type = "ordinal",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalMultiStereotypeLogitRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceOrdinalMultiStereotypeLogitRegr = R6::R6Class("InferenceOrdinalMultiStereotypeLogitRegr",
	inherit = InferenceOrdinalUniStereotypeLogitRegr,
	public = list(

	),

	private = list(
		stereotype_design_matrix = function(){
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

		generate_mod = function(){
			res = fast_stereotype_logit_with_var_cpp(
				X = private$stereotype_design_matrix(),
				y = as.numeric(private$y)
			)
			list(
				b = c(NA, res$b),
				ssq_b_2 = res$ssq_b_1
			)
		}
	)
)
