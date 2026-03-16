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
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignKK14$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "incidence",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 0, 1, 0, 1, 1, 0))
#' infer <- SeqDesignInferenceIncidUnivKKModifiedPoisson$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
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
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
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
