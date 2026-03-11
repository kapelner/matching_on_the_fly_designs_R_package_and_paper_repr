#' Univariate Wilcox Rank-based Regression Compound Inference for KK Designs
#'
#' @description
#' Fits a robust compound estimator for KK matching-on-the-fly designs using rank-based
#' regression (R-estimation) with only the treatment indicator (no additional covariates).
#' For matched pairs, it uses R-estimation on the within-pair response differences.
#' For reservoir subjects, it uses a standard rank-based linear model with the treatment
#' assignment as the sole predictor. This is the univariate (covariate-free) variant.
#'
#' @export
SeqDesignInferenceAllKKWilcoxRegrUnivIVWC = R6::R6Class("SeqDesignInferenceAllKKWilcoxRegrUnivIVWC",
	inherit = SeqDesignInferenceAbstractKKWilcoxRegrIVWC,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design).
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		#' @description
		#' Returns the estimated treatment effect.
		compute_treatment_estimate = function(){
			super$compute_treatment_estimate()
		},

		#' @description
		#' Computes the MLE-based confidence interval.
		#' @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_mle_confidence_interval = function(alpha = 0.05){
			super$compute_mle_confidence_interval(alpha = alpha)
		},

		#' @description
		#' Computes the MLE-based p-value.
		#' @param delta The null difference to test against. For any treatment effect at all this is set to zero (the default).
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			super$compute_mle_two_sided_pval_for_treatment_effect(delta = delta)
		},

		#' @description
		#' Computes the bootstrap confidence interval.
		#' @param alpha The confidence level. Default is 0.05.
		#' @param ... Additional arguments passed to the superclass method.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, ...){
			super$compute_bootstrap_confidence_interval(alpha = alpha, ...)
		}
	),
	private = list(
		include_covariates = function() FALSE
	)
)
