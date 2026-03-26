#' Multivariate Weibull Frailty Combined-Likelihood Inference for KK Designs
#'
#' @description
#' Fits a single joint Weibull frailty model (matched pairs + reservoir as singletons)
#' using the treatment indicator and all recorded covariates.
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
#'   response_type = "survival",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(
#'   c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
#'   c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' infer <- DesignInferenceSurvivalMultiKKWeibullFrailtyCombinedLikelihood$
#'   new(seq_des,
#' verbose = FALSE)
#' infer
#'
DesignInferenceSurvivalMultiKKWeibullFrailtyCombinedLikelihood = R6::R6Class("DesignInferenceSurvivalMultiKKWeibullFrailtyCombinedLikelihood",
	inherit = DesignInferenceAbstractKKWeibullFrailtyCombinedLikelihood,
	public = list(
		#' @description	Initialize the inference object.
		#' @param	des_obj		A SeqDesign object (must be a KK design).
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
		},

		#' @description
		#' Returns the estimated treatment effect.
		compute_treatment_estimate = function(){
			super$compute_treatment_estimate()
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			super$compute_asymp_confidence_interval(alpha = alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta The null difference to test against. For any
		#'   treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			super$compute_asymp_two_sided_pval_for_treatment_effect(delta = delta)
		}
	),
	private = list(
		include_covariates = function() TRUE
	)
)
