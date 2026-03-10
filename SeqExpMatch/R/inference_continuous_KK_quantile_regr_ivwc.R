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
#'   inf = SeqDesignInferenceContinMultKKQuantileRegrIVWC$new(seq_des, tau = 0.75)
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
#' This class requires the \pkg{quantreg} package, which is listed under \code{Suggests}
#' and is not installed automatically with \pkg{SeqExpMatch}. Install it manually with
#' \code{install.packages("quantreg")} before using this class.
#'
#' @export
SeqDesignInferenceContinMultKKQuantileRegrIVWC = R6::R6Class("SeqDesignInferenceContinMultKKQuantileRegrIVWC",
	inherit = SeqDesignInferenceAbstractKKQuantileRegrIVWC,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param	seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param	tau				The quantile level for regression, strictly between 0 and 1. The default \code{tau = 0.5}
		#' 							estimates the median treatment effect. Pass a different value (e.g. \code{tau = 0.25} or
		#' 							\code{tau = 0.75}) to target the corresponding percentile of the treatment effect distribution.
		#' @param	num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs),
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param	verbose			A flag indicating whether messages should be displayed to the user. Default is \code{FALSE}
		initialize = function(seq_des_obj, tau = 0.5, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "continuous")
			super$initialize(seq_des_obj, tau, identity, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		}
	)
)
