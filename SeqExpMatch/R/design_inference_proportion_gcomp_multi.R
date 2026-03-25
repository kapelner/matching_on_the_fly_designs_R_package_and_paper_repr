#' Multivariate G-Computation Mean-Difference Inference for Proportion Responses
#'
#' @description
#' Fits a multivariate fractional-logit working model for a proportion outcome
#' and then estimates the marginal mean difference by standardizing predicted
#' mean proportions under all-treated and all-control assignments over the
#' empirical covariate distribution.
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
#' infer <- DesignInferencePropMultiGCompMeanDiff$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferencePropMultiGCompMeanDiff = R6::R6Class("DesignInferencePropMultiGCompMeanDiff",
	inherit = DesignInferencePropUniGCompMeanDiff,
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
