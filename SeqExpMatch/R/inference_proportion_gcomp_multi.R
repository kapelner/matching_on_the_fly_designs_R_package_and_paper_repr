#' Multivariate G-Computation Mean-Difference Inference for Proportion Responses
#'
#' @description
#' Fits a multivariate fractional-logit working model for a proportion outcome
#' and then estimates the marginal mean difference by standardizing predicted
#' mean proportions under all-treated and all-control assignments over the
#' empirical covariate distribution.
#'
#' @export
SeqDesignInferencePropMultiGCompMeanDiff = R6::R6Class("SeqDesignInferencePropMultiGCompMeanDiff",
	inherit = SeqDesignInferencePropUniGCompMeanDiff,
	public = list(

		#' @description
		#' Initialize the multivariate g-computation mean-difference inference object.
		#' @param seq_des_obj A completed \code{SeqDesign} object with a proportion response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
