#' Univariate Clayton Copula Combined-Likelihood Inference for KK Designs
#'
#' @description
#' Fits a joint copula-based likelihood for KK matching-on-the-fly designs with survival
#' responses using only the treatment indicator. Matched pairs are modeled with a Clayton
#' copula and Weibull AFT margins; reservoir subjects contribute standard Weibull AFT
#' singleton likelihood terms. The treatment effect is reported as a log-time ratio.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignKK14$new(n = nrow(x_dat), response_type = "survival", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(
#'   c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
#'   c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' infer <- SeqDesignInferenceSurvivalUnivKKClaytonCopulaCombinedLikelihood$new(seq_des, verbose
#' = FALSE)
#' infer
#'
SeqDesignInferenceSurvivalUnivKKClaytonCopulaCombinedLikelihood = R6::R6Class("SeqDesignInferenceSurvivalUnivKKClaytonCopulaCombinedLikelihood",
	inherit = SeqDesignInferenceAbstractKKClaytonCopulaCombinedLikelihood,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param seq_des_obj A SeqDesign object (must be a KK design).
		#' @param num_cores Number of CPU cores for parallel processing.
		#' @param verbose Whether to print progress messages.
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
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		compute_mle_confidence_interval = function(alpha = 0.05){
			super$compute_mle_confidence_interval(alpha = alpha)
		},

		#' @description
		#' Computes the MLE-based p-value.
		#' @param delta The null difference to test against. For any
		#'   treatment effect at all this is set to zero (the default).
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			super$compute_mle_two_sided_pval_for_treatment_effect(delta = delta)
		}
	),
	private = list(
		include_covariates = function() FALSE
	)
)
