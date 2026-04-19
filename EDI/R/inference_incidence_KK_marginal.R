#' G-Computation Risk-Difference Inference for KK Designs with Binary Responses
#'
#' Fits an all-subject logistic working model for a KK incidence outcome using
#' treatment and, optionally, all recorded covariates, then estimates the marginal
#' risk difference by standardizing predicted risks under all-treated and
#' all-control assignments over the empirical covariate distribution. Matched
#' pairs are treated as clusters and reservoir subjects are treated as singletons
#' when computing the sandwich covariance.
#'
#' @export
InferenceIncidKKGCompRiskDiff = R6::R6Class("InferenceIncidKKGCompRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceIncidKKGCompAbstract,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			if (private$include_covariates()) {
				private$create_design_matrix()
			} else {
				cbind("(Intercept)" = 1, treatment = private$w)
			}
		},

		get_estimand_type = function() "RD"
	)
)

#' G-Computation Risk-Ratio Inference for KK Designs with Binary Responses
#'
#' Fits a all-subject logistic working model for a KK incidence outcome using
#' treatment and, optionally, all recorded covariates, then estimates the marginal
#' risk ratio by standardizing predicted risks under all-treated and all-control
#' assignments over the empirical covariate distribution. Matched pairs are
#' treated as clusters and reservoir subjects are treated as singletons when
#' computing the sandwich covariance.
#'
#' @export
InferenceIncidKKGCompRiskRatio = R6::R6Class("InferenceIncidKKGCompRiskRatio",
	lock_objects = FALSE,
	inherit = InferenceIncidKKGCompAbstract,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			if (private$include_covariates()) {
				private$create_design_matrix()
			} else {
				cbind("(Intercept)" = 1, treatment = private$w)
			}
		},

		get_estimand_type = function() "RR"
	)
)

#' Modified-Poisson Inference for KK Designs with Binary Responses
#'
#' Fits an all-subject modified-Poisson working model for incidence outcomes under
#' a KK matching-on-the-fly design using treatment and, optionally, all recorded
#' covariates as predictors. Matched pairs are treated as clusters and reservoir
#' subjects are treated as singleton clusters when computing the sandwich
#' covariance.
#'
#' @export
InferenceIncidKKModifiedPoisson = R6::R6Class("InferenceIncidKKModifiedPoisson",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKModifiedPoisson,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			if (private$include_covariates()) {
				private$create_design_matrix()
			} else {
				cbind("(Intercept)" = 1, treatment = private$w)
			}
		}
	)
)
