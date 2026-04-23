#' Zero-Inflated Poisson Regression Inference for Count Responses
#'
#' Fits a zero-inflated Poisson regression for count responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceCountZeroInflatedPoisson = R6::R6Class("InferenceCountZeroInflatedPoisson",
	lock_objects = FALSE,
	inherit = InferenceCountZeroAugmentedPoissonAbstract,
	public = list(
	),
	private = list(
		za_family = function() stats::poisson(link = "log"),
		za_description = function() "Zero-Inflated Poisson"
	)
)

#' Zero-Inflated Negative Binomial Regression Inference for Count Responses
#'
#' Fits a zero-inflated negative binomial regression for count responses using
#' the treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceCountZeroInflatedNegBin = R6::R6Class("InferenceCountZeroInflatedNegBin",
	lock_objects = FALSE,
	inherit = InferenceCountZeroAugmentedPoissonAbstract,
	public = list(
		#' @description
		#' Initialize a zero-inflated negative binomial inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp
		#'   implementation. If \code{FALSE}, use \pkg{glmmTMB}.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg Optimization algorithm. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, optimization_alg = NULL){
			super$initialize(des_obj, model_formula = model_formula, use_rcpp = use_rcpp, verbose = verbose, optimization_alg = optimization_alg)
		}
	),
	private = list(
		za_family = function() glmmTMB::nbinom2(link = "log"),
		za_description = function() "Zero-Inflated Negative Binomial"
	)
)
