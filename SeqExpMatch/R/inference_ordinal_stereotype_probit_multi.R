#' Multivariate Stereotype Probit Inference for Ordinal Responses
#
#' @description
#' Stereotype probit inference for ordinal responses using treatment plus observed covariates.
#'
#' @inherit InferenceRand methods
#' @inherit InferenceBoot methods
#' @inherit InferenceAsymp methods
#' @inherit InferenceRandCI methods
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
#' infer <- InferenceOrdinalMultiStereotypeProbitRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceOrdinalMultiStereotypeProbitRegr = R6::R6Class("InferenceOrdinalMultiStereotypeProbitRegr",
	inherit = InferenceOrdinalUniStereotypeProbitRegr,
	public = list(

		#' @description
		#' Initialize a multivariate stereotype probit inference object for a
		#' completed ordinal sequential design.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with an ordinal response.
		#' @param num_cores Number of CPU cores for bootstrap/randomization helpers.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
		}
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
		}
	)
)
