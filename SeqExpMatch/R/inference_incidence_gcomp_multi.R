#' Multivariate G-Computation Risk-Difference Inference for Binary Responses
#'
#' @description
#' Fits a logistic working model for an incidence outcome using treatment and all
#' recorded covariates, then estimates the marginal risk difference by
#' standardizing predicted risks under all-treated and all-control assignments
#' over the empirical covariate distribution.
#'
#' @export
SeqDesignInferenceIncidMultiGCompRiskDiff = R6::R6Class("SeqDesignInferenceIncidMultiGCompRiskDiff",
	inherit = SeqDesignInferenceIncidUnivGCompRiskDiff,
	public = list(

		#' @description
		#' Initialize the multivariate g-computation RD inference object.
		#' @param seq_des_obj A completed \code{SeqDesign} object with an incidence response.
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

#' Multivariate G-Computation Risk-Ratio Inference for Binary Responses
#'
#' @description
#' Fits a logistic working model for an incidence outcome using treatment and all
#' recorded covariates, then estimates the marginal risk ratio by standardizing
#' predicted risks under all-treated and all-control assignments over the
#' empirical covariate distribution.
#'
#' @details
#' The point estimate is returned on the risk-ratio scale. Confidence intervals
#' and p-values use the delta method on the log-risk-ratio scale and then map
#' back to the risk-ratio scale.
#'
#' @export
SeqDesignInferenceIncidMultiGCompRiskRatio = R6::R6Class("SeqDesignInferenceIncidMultiGCompRiskRatio",
	inherit = SeqDesignInferenceIncidUnivGCompRiskRatio,
	public = list(

		#' @description
		#' Initialize the multivariate g-computation RR inference object.
		#' @param seq_des_obj A completed \code{SeqDesign} object with an incidence response.
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
