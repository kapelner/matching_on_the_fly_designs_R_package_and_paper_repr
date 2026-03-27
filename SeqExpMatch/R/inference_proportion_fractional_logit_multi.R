#' Multivariate Fractional Logit Inference for Proportion Responses
#'
#' @description
#' Fits a fractional logit model for proportion responses using the treatment
#' indicator and all recorded covariates, with sandwich-robust variance. The
#' treatment effect is reported on the log-odds scale.
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
#'   response_type = "proportion",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0.10, 0.25, 0.20, 0.40, 0.35, 0.55, 0.60, 0.75))
#' infer <- InferencePropMultiFractionalLogit$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferencePropMultiFractionalLogit = R6::R6Class("InferencePropMultiFractionalLogit",
	inherit = InferencePropUniFractionalLogit,
	public = list(

		#' @description
		#' Initialize a multivariate fractional-logit inference object for a completed
		#' design with a proportion response.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with a proportion response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
		}
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
