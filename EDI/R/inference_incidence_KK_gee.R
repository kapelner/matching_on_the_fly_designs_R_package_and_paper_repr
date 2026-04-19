#' GEE Inference for KK Designs with Binary Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for binary (incidence) responses under a KK matching-on-the-fly design using
#' the treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @details
#' This class requires the \pkg{geepack} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{geepack} before using this class.
#'
#' @export
InferenceIncidKKGEE = R6::R6Class("InferenceIncidKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(
	),
	private = list(
		gee_response_type = function() "incidence",
		gee_family        = function() stats::binomial(link = "logit")
	)
)
