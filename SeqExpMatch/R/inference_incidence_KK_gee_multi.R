#' Multivariate GEE Inference for KK Designs with Binary Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for binary (incidence) responses under a KK matching-on-the-fly design using the
#' treatment indicator and all recorded covariates as predictors. Matched pairs are
#' treated as clusters (with exchangeable correlation structure); reservoir subjects
#' each form their own singleton cluster. Unlike conditional-logit-only matched-pair
#' analyses, all subjects (matched and reservoir) are included. Inference is based on
#' sandwich-robust standard errors, so the test statistic is Z-distributed.
#'
#' @details
#' This class requires the \pkg{geepack} package, which is listed in Suggests
#' and is not installed automatically with \pkg{SeqExpMatch}.
#' Install \pkg{geepack} before using this class.
#'
#' @export
InferenceIncidMultiKKGEE = R6::R6Class("InferenceIncidMultiKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(




	),
	private = list(
		gee_response_type = function() "incidence",
		gee_family        = function() binomial(link = "logit")
	)
)
