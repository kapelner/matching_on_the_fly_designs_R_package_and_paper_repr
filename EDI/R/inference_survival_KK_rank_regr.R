#' Survival Rank-based Regression (AFT) Compound Inference for KK Designs
#'
#' Fits a robust compound estimator for KK matching-on-the-fly designs with survival
#' responses using rank-based estimating equations.
#'
#' @export
InferenceSurvivalKKRankRegrIVWC = R6::R6Class("InferenceSurvivalKKRankRegrIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKSurvivalRankRegrIVWC,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			super$initialize(des_obj, include_covariates, verbose)
		}
	)
)
