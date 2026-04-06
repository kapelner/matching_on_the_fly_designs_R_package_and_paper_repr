#' Univariate GEE Inference for KK Designs with Proportion Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for proportion (continuous values in (0, 1)) responses under a KK
#' matching-on-the-fly design using only the treatment indicator as a predictor
#' (intercept + treatment). A binomial logit working model is used; although
#' designed for binary data, \code{geepack} accepts continuous (0, 1) responses
#' with this family and it is more numerically stable than \code{quasibinomial}
#' (which lacks a proper log-likelihood and can cause QIC computation failures in
#' \code{geepack}). The sandwich-robust standard errors are valid regardless of
#' working-model misspecification and already absorb any overdispersion. Matched
#' pairs are treated as clusters (with exchangeable correlation); reservoir subjects
#' each form their own singleton cluster. The treatment estimate is the log-odds
#' ratio of the mean proportion; inference uses Wald Z-statistics based on the
#' robust SE.
#'
#' @details
#' This class requires the \pkg{geepack} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{geepack} before using this class.
#'
#' @export
InferencePropUnivKKGEE = R6::R6Class("InferencePropUnivKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(




	),
	private = list(
		gee_response_type = function() "proportion",
		gee_family        = function() binomial(link = "logit"),
		gee_predictors_df = function() data.frame(w = private$w)
	)
)
