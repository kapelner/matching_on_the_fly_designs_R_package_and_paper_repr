#' Multivariate Modified-Poisson Inference for KK Designs with Binary Responses
#'
#' @description
#' Fits an all-subject modified-Poisson working model for incidence outcomes under
#' a KK matching-on-the-fly design using treatment and all recorded covariates as
#' predictors. Matched pairs are treated as clusters and reservoir subjects are
#' treated as singleton clusters when computing the sandwich covariance, so the
#' estimated treatment effect is a log risk ratio with cluster-robust inference.
#'
#' @export
SeqDesignInferenceIncidMultiKKModifiedPoisson = R6::R6Class("SeqDesignInferenceIncidMultiKKModifiedPoisson",
	inherit = SeqDesignInferenceIncidUnivKKModifiedPoisson,
	public = list(

		#' @description
		#' Initialize a multivariate modified-Poisson inference object for a completed
		#' KK design with a binary response.
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
