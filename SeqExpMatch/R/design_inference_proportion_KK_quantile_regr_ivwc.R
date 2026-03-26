#' Quantile Regression Compound Estimator for KK Matching-on-the-Fly Designs (Proportion
#' Outcomes)
#'
#' @description
#' A variance-weighted compound quantile regression estimator for KK matching-on-the-fly
#' designs with proportion responses. Inference is performed on the \strong{logit (log-odds)
#' scale}: responses \eqn{y \in (0,1)} are transformed via \eqn{\text{logit}(y) = \log(y/(1-y))}
#' before quantile regression.
#'
#' The estimator combines:
#' \enumerate{
#'   \item Quantile regression on logit-scale within-pair differences
#'     \eqn{\text{logit}(y_T) - \text{logit}(y_C)} (matched pairs)
#'   \item Quantile regression of \eqn{\text{logit}(y)} on treatment and covariates (reservoir)
#' }
#' using the same variance-weighted combination logic as the OLS compound estimator.
#'
#' The estimated treatment effect is a \strong{log-odds-ratio shift} at quantile \code{tau}.
#' At \code{beta_T = 1} (one log-odds-ratio unit of treatment effect), the population
#' treatment effect on the logit scale is exactly 1, so no \code{skip_ci} is needed.
#'
#' \strong{Default quantile: \code{tau = 0.5} (median regression).}
#' To target a different quantile — for example the 25th or 75th percentile — pass
#' \code{tau = 0.25} or \code{tau = 0.75} to the constructor:
#' \preformatted{
#'   inf = DesignInferencePropMultiKKQuantileRegrIVWC$
#'   new(seq_des, tau = 0.75)
#' }
#' Any value strictly between 0 and 1 is accepted.
#'
#' Standard errors use Powell's "nid" sandwich estimator (non-iid), falling back to "iid"
#' on failure. Asymptotic z-based inference is used throughout.
#'
#' This class requires the \pkg{quantreg} package, which is listed in Suggests
#' and is not installed automatically with \pkg{SeqExpMatch}.
#' Install \pkg{quantreg} before using this class.
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
#'   response_type = "proportion",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0.10, 0.25, 0.20, 0.40, 0.35, 0.55, 0.60, 0.75))
#' infer <- DesignInferencePropMultiKKQuantileRegrIVWC$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferencePropMultiKKQuantileRegrIVWC = R6::R6Class("DesignInferencePropMultiKKQuantileRegrIVWC",
	inherit = DesignInferenceAbstractKKQuantileRegrIVWC,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj A SeqDesign object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param tau                             The quantile level for regression on the logit
		#'   scale, strictly between 0 and 1.
		#' 							The default \code{tau = 0.5} estimates the median log-odds-ratio treatment effect.
		#' 							Pass a different value (e.g. \code{tau = 0.25} or \code{tau = 0.75}) to target a
		#' 							different percentile of the treatment effect distribution.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#' 							and bootstrap resampling. The default is 1 for serial computation.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{FALSE}.
		initialize = function(des_obj, tau = 0.5, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "proportion")
			super$initialize(des_obj, tau, qlogis, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
			assertNumeric(private$y, any.missing = FALSE, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
		},

		#' @description
		#' Computes the randomization-based confidence interval via Zhang's combined test.
		#' @param alpha The confidence level is 1 - \code{alpha}.
		#' @param r Number of random sign-flips / permutations.
		#' @param pval_epsilon Bisection convergence tolerance.
		#' @param show_progress Ignored.
		compute_confidence_interval_rand = function(alpha = 0.05, r = 499, pval_epsilon = 0.005, show_progress = TRUE){
			super$compute_confidence_interval_rand(
				alpha = alpha,
				r = r,
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
