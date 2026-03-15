#' Multivariate Risk-Difference Regression for Binary Responses
#'
#' @description
#' Fits a direct risk-difference estimator for binary (incidence) responses
#' under non-KK designs using treatment and all recorded covariates in a linear
#' probability model with HC2 heteroskedasticity-robust variance. The treatment
#' effect is reported on the risk-difference scale.
#'
#' @export
SeqDesignInferenceIncidMultiRiskDiff = R6::R6Class("SeqDesignInferenceIncidMultiRiskDiff",
	inherit = SeqDesignInferenceIncidUnivRiskDiff,
	public = list(

		#' @description
		#' Initialize a multivariate risk-difference inference object for a
		#' completed non-KK design with a binary response.
		#' @param seq_des_obj A completed non-KK \code{SeqDesign} object with an
		#'   incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and
		#'   randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		#' @description
		#' Computes the direct risk-difference estimate of the treatment effect.
		compute_treatment_estimate = function(){
			super$compute_treatment_estimate()
		}
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
