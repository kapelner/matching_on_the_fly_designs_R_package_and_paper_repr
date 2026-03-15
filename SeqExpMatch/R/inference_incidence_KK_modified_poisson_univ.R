#' Univariate Modified-Poisson Inference for KK Designs with Binary Responses
#'
#' @description
#' Fits an all-subject modified-Poisson working model for incidence outcomes under
#' a KK matching-on-the-fly design using only the treatment indicator as a
#' predictor. Matched pairs are treated as clusters and reservoir subjects are
#' treated as singleton clusters when computing the sandwich covariance, so the
#' estimated treatment effect is a log risk ratio with cluster-robust inference.
#'
#' @export
SeqDesignInferenceIncidUnivKKModifiedPoisson = R6::R6Class("SeqDesignInferenceIncidUnivKKModifiedPoisson",
	inherit = SeqDesignInferenceAbstractKKModifiedPoisson,
	public = list(

		#' @description
		#' Initialize a univariate modified-Poisson inference object for a completed
		#' KK design with a binary response.
		#' @param seq_des_obj A completed KK \code{SeqDesign} object with an incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		#' @description
		#' Computes the modified-Poisson estimate of the treatment effect on the
		#' log-risk-ratio scale.
		compute_treatment_estimate = function(){
			super$compute_treatment_estimate()
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval using the cluster-robust
		#' sandwich standard error.
		#' @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_mle_confidence_interval = function(alpha = 0.05){
			super$compute_mle_confidence_interval(alpha = alpha)
		},

		#' @description
		#' Computes a two-sided p-value for the treatment effect.
		#' @param delta The null treatment effect on the log-risk-ratio scale.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			super$compute_mle_two_sided_pval_for_treatment_effect(delta = delta)
		}
	),

	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		}
	)
)
