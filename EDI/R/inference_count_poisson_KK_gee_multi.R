#' Multivariate GEE Inference for KK Designs with Count Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for count responses under a KK matching-on-the-fly design using the treatment
#' indicator and all recorded covariates as predictors. A Poisson working model with
#' a log link is used; the sandwich-robust standard errors remain valid even under
#' overdispersion. Matched pairs are treated as clusters (with exchangeable
#' correlation); reservoir subjects each form their own singleton cluster. The
#' treatment estimate is the log incidence-rate ratio (log-IRR).
#'
#' @details
#' This class requires the \pkg{geepack} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{geepack} before using this class.
#'
#' @export
InferenceCountPoissonMultiKKGEE = R6::R6Class("InferenceCountPoissonMultiKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(




	),
	private = list(
		gee_response_type = function() "count",
		gee_family        = function() poisson(link = "log")
	)
)
