#' Quantile Regression Combined-Likelihood Compound Estimator for KK Designs (Proportion)
#'
#' @description
#' Fits the combined stacked quantile regression (matched-pair differences + reservoir)
#' using the treatment indicator and all recorded covariates for proportion responses.
#' Responses are transformed via logit before regression; the estimated treatment
#' effect is a log-odds-ratio shift at quantile \code{tau}.
#' Minimises the joint check-function loss over both data sources simultaneously.
#' Inference is based on the stacked combined-likelihood quantile-regression fit.
#'
#' @export
SeqDesignInferencePropMultiKKQuantileRegrCombinedLikelihood = R6::R6Class("SeqDesignInferencePropMultiKKQuantileRegrCombinedLikelihood",
	inherit = SeqDesignInferenceAbstractKKQuantileRegrCombinedLikelihood,
	public = list(
		#' @description	Initialize the inference object.
		#' @param	seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param	tau				The quantile level on the logit scale, strictly between 0 and 1. Default is 0.5.
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, tau = 0.5, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "proportion")
			super$initialize(seq_des_obj, tau, qlogis, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
			assertNumeric(private$y, any.missing = FALSE, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
		},

		#' @description
		#' Computes the randomization-based confidence interval via Zhang's combined test.
		#' @param alpha The confidence level is 1 - \code{alpha}.
		#' @param nsim_exact_test Number of random sign-flips / permutations.
		#' @param pval_epsilon Bisection convergence tolerance.
		#' @param show_progress Ignored.
		compute_confidence_interval_rand = function(alpha = 0.05, nsim_exact_test = 499, pval_epsilon = 0.005, show_progress = TRUE){
			super$compute_confidence_interval_rand(
				alpha = alpha,
				nsim_exact_test = nsim_exact_test,
				pval_epsilon = pval_epsilon,
				show_progress = show_progress
			)
		},

		#' @description
		#' Returns the estimated treatment effect.
		compute_treatment_estimate = function() super$compute_treatment_estimate(),

		#' @description
		#' Computes the MLE-based confidence interval.
		#' @param alpha Significance level; default 0.05 gives a 95 percent CI.
		compute_mle_confidence_interval = function(alpha = 0.05){
			super$compute_mle_confidence_interval(alpha = alpha)
		},

		#' @description
		#' Computes the MLE-based p-value.
		#' @param delta Null value; default 0.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			super$compute_mle_two_sided_pval_for_treatment_effect(delta = delta)
		}
	)
)
