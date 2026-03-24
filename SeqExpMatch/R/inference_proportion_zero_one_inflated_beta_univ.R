#' Univariate Zero/One-Inflated Beta Inference for Proportion Responses
#'
#' @description
#' Fits a zero/one-inflated beta regression for proportion responses using only the
#' treatment indicator in the beta mean submodel. Exact 0 and exact 1 values are
#' modeled with separate inflation masses, and interior values are modeled with a
#' beta distribution. The reported treatment effect is the treatment coefficient
#' from the beta mean submodel, on the logit scale.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignBernoulli$
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
#' infer <- SeqDesignInferencePropUniZeroOneInflatedBetaRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
SeqDesignInferencePropUniZeroOneInflatedBetaRegr = R6::R6Class("SeqDesignInferencePropUniZeroOneInflatedBetaRegr",
	inherit = SeqDesignInferencePropZeroOneInflatedBetaAbstract,
	public = list(

		#' @description
		#' Initialize a zero/one-inflated beta inference object for a completed design
		#' with a proportion response.
		#' @param seq_des_obj A completed \code{SeqDesign} object with a proportion response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	)
)
