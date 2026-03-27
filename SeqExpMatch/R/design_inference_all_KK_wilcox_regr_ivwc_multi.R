#' Multivariate Wilcox Rank-based Regression Compound Inference for KK Designs
#'
#' @description
#' Fits a robust compound estimator for KK matching-on-the-fly designs using rank-based
#' regression (R-estimation) with the treatment indicator and all recorded covariates.
#' For matched pairs, it uses R-estimation on the within-pair differences of responses
#' and covariates. For reservoir subjects, it uses a standard rank-based linear model
#' including both treatment and covariates. This is the multivariate (covariate-adjusted)
#' variant.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
DesignInferenceAllKKWilcoxRegrMultiIVWC = R6::R6Class("DesignInferenceAllKKWilcoxRegrMultiIVWC",
	inherit = DesignInferenceAbstractKKWilcoxRegrIVWC,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param	des_obj		A SeqDesign object (must be a KK design).
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		#' @examples
		#' set.seed(1)
		#' x_dat <- data.frame(
		#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
		#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
		#' )
		#' seq_des <- SeqDesignKK14$new(n = nrow(x_dat), response_type = "continuous", verbose =
		#' FALSE)
		#' for (i in seq_len(nrow(x_dat))) {
		#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
		#' }
		#' seq_des$add_all_subject_responses(c(1.2, 0.9, 1.5, 1.8, 2.1, 1.7, 2.6, 2.2))
		#' infer <- DesignInferenceAllKKWilcoxRegrMultiIVWC$new(seq_des, verbose = FALSE)
		#' infer
		#'
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
		include_covariates = function() TRUE
	)
)
