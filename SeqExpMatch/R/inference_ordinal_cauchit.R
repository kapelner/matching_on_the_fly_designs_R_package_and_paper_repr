#' Cumulative Cauchit Inference for Ordinal Responses
#'
#' @description
#' Cumulative Cauchit model inference for ordinal responses.
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
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalUniCauchitRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalUniCauchitRegr = R6::R6Class("InferenceOrdinalUniCauchitRegr",
	inherit = InferenceMLEorKMforGLMs,
	public = list(
		#' @description
		#' Initialize a sequential experimental design estimation and test object.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects are assigned and
		#'   response y is recorded within.
		#' @param num_cores The number of CPU cores to use.
		#' @param verbose A flag indicating whether messages should be displayed.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the appropriate estimate
		#' @return The setting-appropriate numeric estimate of the treatment effect
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		}
	),

	private = list(
		generate_mod = function(){
			Xmm = matrix(private$w, ncol = 1)
			colnames(Xmm) = c("treatment")
			res = fast_ordinal_cauchit_regression_with_var_cpp(X = Xmm, y = as.numeric(private$y))
			list(
				b = c(NA, res$b[1]),
				ssq_b_2 = res$ssq_b_2
			)
		}
	)
)

#' Multivariate Cumulative Cauchit Inference for Ordinal Responses
#'
#' @description
#' Cumulative Cauchit model inference for ordinal responses with
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
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalMultiCauchitRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalMultiCauchitRegr = R6::R6Class("InferenceOrdinalMultiCauchitRegr",
	inherit = InferenceOrdinalUniCauchitRegr,
	public = list(
	),

	private = list(
		cauchit_design_matrix = function(){
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
			res = fast_ordinal_cauchit_regression_with_var_cpp(
				X = private$cauchit_design_matrix(),
				y = as.numeric(private$y)
			)
			list(
				b = c(NA, res$b),
				ssq_b_2 = res$ssq_b_2
			)
		}
	)
)
