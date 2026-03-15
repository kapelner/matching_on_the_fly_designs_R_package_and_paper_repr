#' Multivariate Fractional Logit Inference for Proportion Responses
#'
#' @description
#' Fits a fractional logit model for proportion responses using the treatment
#' indicator and all recorded covariates, with sandwich-robust variance. The
#' treatment effect is reported on the log-odds scale.
#'
#' @export
SeqDesignInferencePropMultiFractionalLogit = R6::R6Class("SeqDesignInferencePropMultiFractionalLogit",
	inherit = SeqDesignInferencePropUniFractionalLogit,
	public = list(

		#' @description
		#' Initialize a multivariate fractional-logit inference object for a completed
		#' design with a proportion response.
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
