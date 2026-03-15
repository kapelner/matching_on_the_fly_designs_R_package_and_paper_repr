#' Multivariate G-Computation Risk-Difference Inference for KK Designs with Binary Responses
#'
#' @description
#' Fits an all-subject logistic working model for a KK incidence outcome using
#' treatment and all recorded covariates, then estimates the marginal risk
#' difference by standardizing predicted risks under all-treated and all-control
#' assignments over the empirical covariate distribution. Matched pairs are
#' treated as clusters and reservoir subjects are treated as singletons when
#' computing the sandwich covariance.
#'
#' @export
SeqDesignInferenceIncidMultiKKGCompRiskDiff = R6::R6Class("SeqDesignInferenceIncidMultiKKGCompRiskDiff",
	inherit = SeqDesignInferenceIncidUnivKKGCompRiskDiff,
	public = list(
		#' @description
		#' Initialize the multivariate KK g-computation RD inference object.
		#' @param seq_des_obj A completed KK \code{SeqDesign} object with an incidence response.
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

#' Multivariate G-Computation Risk-Ratio Inference for KK Designs with Binary Responses
#'
#' @description
#' Fits an all-subject logistic working model for a KK incidence outcome using
#' treatment and all recorded covariates, then estimates the marginal risk ratio
#' by standardizing predicted risks under all-treated and all-control assignments
#' over the empirical covariate distribution. Matched pairs are treated as
#' clusters and reservoir subjects are treated as singletons when computing the
#' sandwich covariance.
#'
#' @details
#' The point estimate is returned on the risk-ratio scale. Confidence intervals
#' and p-values use the delta method on the log-risk-ratio scale and then map
#' back to the risk-ratio scale.
#'
#' @export
SeqDesignInferenceIncidMultiKKGCompRiskRatio = R6::R6Class("SeqDesignInferenceIncidMultiKKGCompRiskRatio",
	inherit = SeqDesignInferenceIncidUnivKKGCompRiskRatio,
	public = list(
		#' @description
		#' Initialize the multivariate KK g-computation RR inference object.
		#' @param seq_des_obj A completed KK \code{SeqDesign} object with an incidence response.
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
