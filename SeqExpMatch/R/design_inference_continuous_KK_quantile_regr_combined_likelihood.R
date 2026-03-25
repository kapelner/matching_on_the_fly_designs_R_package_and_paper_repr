#' Quantile Regression Combined-Likelihood Compound Estimator for KK Designs (Continuous)
#'
#' @description
#' Fits the combined stacked quantile regression (matched-pair differences + reservoir)
#' using the treatment indicator and all recorded covariates for continuous responses.
#' Minimises the joint check-function loss over both data sources simultaneously.
#' Inference is based on the stacked combined-likelihood quantile-regression fit.
#'
#' @export
DesignInferenceContinMultKKQuantileRegrCombinedLikelihood = R6::R6Class("DesignInferenceContinMultKKQuantileRegrCombinedLikelihood",
	inherit = DesignInferenceAbstractKKQuantileRegrCombinedLikelihood,
	public = list(
		#' @description	Initialize the inference object.
		#' @param seq_des_obj A SeqDesign object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param	tau				The quantile level for regression, strictly between 0 and 1. Default is 0.5.
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
		#' infer <- DesignInferenceContinMultKKQuantileRegrCombinedLikelihood$new(seq_des, verbose
		#' = FALSE)
		#' infer
		#'
		initialize = function(seq_des_obj, tau = 0.5, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "continuous")
			super$initialize(seq_des_obj, tau, identity, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
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
		#' Computes the asymptotic confidence interval.
		#' @param alpha Significance level; default 0.05 gives a 95 percent CI.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			super$compute_asymp_confidence_interval(alpha = alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta Null value; default 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			super$compute_asymp_two_sided_pval_for_treatment_effect(delta = delta)
		}
	)
)
