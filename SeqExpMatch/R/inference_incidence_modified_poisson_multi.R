#' Multivariate Modified Poisson Inference for Binary Responses
#'
#' @description
#' Fits the classic modified Poisson model for binary (incidence) responses under
#' non-KK designs using treatment and all recorded covariates as predictors in a
#' Poisson log-link working model with Huber-White sandwich variance. The
#' treatment effect is reported on the log-risk-ratio scale.
#'
#' @export
SeqDesignInferenceIncidMultiModifiedPoisson = R6::R6Class("SeqDesignInferenceIncidMultiModifiedPoisson",
	inherit = SeqDesignInferenceIncidUnivModifiedPoisson,
	public = list(

		#' @description
		#' Initialize a multivariate modified-Poisson inference object for a
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
		#' Computes the modified-Poisson estimate of the treatment effect on the
		#' log-risk-ratio scale.
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
