#' Quantile Regression Compound Estimator for KK Matching-on-the-Fly Designs
#'
#' @description
#' A variance-weighted compound quantile regression estimator for KK matching-on-the-fly
#' designs with continuous responses. The estimator combines:
#' \enumerate{
#'   \item Quantile regression on within-pair differences (matched pairs)
#'   \item Quantile regression on reservoir subjects (treatment vs control)
#' }
#' using the same variance-weighted combination logic as the OLS compound estimator.
#'
#' \strong{Default quantile: \code{tau = 0.5} (median regression).}
#' At \code{tau = 0.5} this estimates the median treatment effect, which is the canonical
#' nonparametric location estimator and is more robust to outliers and heavy-tailed
#' response distributions than the OLS mean-based estimator. To target a different
#' quantile of the treatment effect distribution — for example the 25th or 75th
#' percentile — pass \code{tau = 0.25} or \code{tau = 0.75} to the constructor:
#' \preformatted{
#'   inf = DesignInferenceContinMultKKQuantileRegrIVWC$
#'   new(seq_des, tau = 0.75)
#' }
#' Any value strictly between 0 and 1 is accepted.
#'
#' Standard errors use Powell's "nid" sandwich estimator (non-iid), which is more robust
#' than the "iid" (constant-density) assumption; the implementation falls back to "iid"
#' on failure. Asymptotic z-based inference is used throughout.
#'
#' The randomization-based confidence interval is inherited from the base class and is
#' valid for location-shift models at all quantiles: shifting y by delta maps the
#' tau-th quantile treatment effect to delta under the null.
#'
#' This class requires the \pkg{quantreg} package, which is listed in Suggests
#' and is not installed automatically with \pkg{SeqExpMatch}.
#' Install \pkg{quantreg} before using this class.
#'
#' @export
DesignInferenceContinMultKKQuantileRegrIVWC = R6::R6Class("DesignInferenceContinMultKKQuantileRegrIVWC",
	inherit = DesignInferenceAbstractKKQuantileRegrIVWC,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param seq_des_obj A SeqDesign object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param tau                             The quantile level for regression, strictly between
		#'   0 and 1. The default \code{tau = 0.5}
		#'                                                         estimates the median treatment
		#' effect. Pass a different value (e.g. \code{tau = 0.25} or
		#'                                                         \code{tau = 0.75}) to target the
		#' corresponding percentile of the treatment effect distribution.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#'   The default is 1 for serial computation. For simple
		#'   estimators (e.g. mean difference and KK compound),
		#'   parallelization is achieved with zero-overhead C++ OpenMP.
		#'   For complex models (e.g. GLMs),
		#'   parallelization falls back to R's
		#'   \code{parallel::mclapply}, which incurs
		#'   session-forking overhead.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{FALSE}.
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
		#' infer <- DesignInferenceContinMultKKQuantileRegrIVWC$new(seq_des, verbose = FALSE)
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
		#' @param delta The null difference to test against. Default is zero.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			super$compute_asymp_two_sided_pval_for_treatment_effect(delta = delta)
		}
	)
)
