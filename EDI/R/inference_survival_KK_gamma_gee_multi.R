#' Multivariate GEE Inference for KK Designs with Survival Response (no censoring)
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for uncensored survival (time-to-event) responses under a KK matching-on-the-fly
#' design using the treatment indicator and all recorded covariates as predictors.
#' A Gamma working model with a log link is used, which is appropriate for positive
#' continuous survival times. The sandwich-robust standard errors from GEE are valid
#' regardless of Gamma misspecification. Matched pairs are treated as clusters (with
#' exchangeable correlation); reservoir subjects each form their own singleton cluster.
#' The treatment estimate is the log ratio of mean survival times (log-MTR); inference
#' uses Wald Z-statistics based on the robust SE.
#'
#' @details
#' Censored observations are not supported; an error is raised at initialization if
#' any censoring is detected. This class requires the \pkg{geepack} package, which is
#' listed under \code{Suggests} and is not installed automatically with
#' \pkg{EDI}. Install \pkg{geepack} manually
#' before using this class.
#'
#' @export
InferenceSurvivalMultiKKGammaGEE = R6::R6Class("InferenceSurvivalMultiKKGammaGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(
		#' @description
		#' Initialize the Inference object.
		#'
		#' @param des_obj The design object.
		#' @param verbose If TRUE, print additional information.
		initialize = function(des_obj, verbose = FALSE) {
			assertResponseType(des_obj$get_response_type(), "survival")
			super$initialize(des_obj, verbose)
		}




	),
	private = list(
		gee_response_type = function() "survival",
		gee_family        = function() Gamma(link = "log")
	)
)
