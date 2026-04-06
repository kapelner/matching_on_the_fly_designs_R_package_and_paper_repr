#' Univariate GEE Inference for KK Designs with Survival Response (no censoring)
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for uncensored survival (time-to-event) responses under a KK matching-on-the-fly
#' design using only the treatment indicator as a predictor (intercept + treatment).
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
InferenceSurvivalUnivKKGEE = R6::R6Class("InferenceSurvivalUnivKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(




	),
	private = list(
		gee_response_type = function() "survival",
		gee_family        = function() Gamma(link = "log"),
		gee_predictors_df = function() data.frame(w = private$w)
	)
)
