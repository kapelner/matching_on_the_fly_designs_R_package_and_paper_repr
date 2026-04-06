#' Quantile Regression Compound Estimator for KK Matching-on-the-Fly Designs (Proportion
#' Outcomes)
#'
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
#'   inf = InferencePropMultiKKQuantileRegrIVWC$
#'   new(seq_des, tau = 0.75)
#' }
#' Any value strictly between 0 and 1 is accepted.
#'
#' Standard errors use Powell's "nid" sandwich estimator (non-iid), falling back to "iid"
#' on failure. Asymptotic z-based inference is used throughout.
#'
#' This class requires the \pkg{quantreg} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{quantreg} before using this class.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneKK14$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "proportion",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0.10, 0.25, 0.20, 0.40, 0.35, 0.55, 0.60, 0.75))
#' infer <- InferencePropMultiKKQuantileRegrIVWC$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferencePropMultiKKQuantileRegrIVWC = R6::R6Class("InferencePropMultiKKQuantileRegrIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKQuantileRegrIVWC,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param tau                             The quantile level for regression on the logit
		#'   scale, strictly between 0 and 1.
		#' 							The default \code{tau = 0.5} estimates the median log-odds-ratio treatment effect.
		#' 							Pass a different value (e.g. \code{tau = 0.25} or \code{tau = 0.75}) to target a
		#' 							different percentile of the treatment effect distribution.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{FALSE}.
		initialize = function(des_obj, tau = 0.5,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "proportion")
			super$initialize(des_obj, tau, qlogis, verbose)
			assertNoCensoring(private$any_censoring)
			private$y = .sanitize_proportion_response(private$y, interior = TRUE)
			assertNumeric(private$y, any.missing = FALSE, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			# Rebuild KK summary data after sanitizing the response, otherwise the
			# superclass cache can retain raw 0/1 values and qlogis() will produce
			# non-finite values during the quantile fit and its randomization refits.
			private$cached_values$KKstats = NULL
			private$compute_basic_match_data()
		}




	)
)
