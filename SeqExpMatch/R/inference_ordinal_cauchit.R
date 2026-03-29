#' Cumulative Cauchit Inference for Ordinal Responses
#'
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
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(
		#' @description
		#' Initialize a sequential experimental design estimation and test object.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects are assigned and
		#'   response y is recorded within.
		#' @param num_cores The number of CPU cores to use.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param make_fork_cluster Whether to use a fork cluster for parallelization.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
		}
	),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			Xmm = matrix(private$w, ncol = 1)
			colnames(Xmm) = c("treatment")
			
			if (estimate_only) {
				res = fast_ordinal_cauchit_regression_cpp(X = Xmm, y = as.numeric(private$y))
				return(list(
					b = c(NA, res$b[1]),
					ssq_b_2 = NA_real_
				))
			}

			res = fast_ordinal_cauchit_regression_with_var_cpp(X = Xmm, y = as.numeric(private$y))
			list(
				b = c(NA, res$b[1]),
				ssq_b_2 = res$ssq_b_2
			)
		}
	)
)
InferenceOrdinalMultiCauchitRegr = R6::R6Class("InferenceOrdinalMultiCauchitRegr",
	lock_objects = FALSE,
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

		generate_mod = function(estimate_only = FALSE){
			if (estimate_only) {
				res = fast_ordinal_cauchit_regression_cpp(
					X = private$cauchit_design_matrix(),
					y = as.numeric(private$y)
				)
				return(list(
					b = c(NA, res$b),
					ssq_b_2 = NA_real_
				))
			}

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
