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
