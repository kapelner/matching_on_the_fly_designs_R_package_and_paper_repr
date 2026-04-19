#' Stratified Cox / Standard Cox Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with survival responses.
#' For matched pairs, it uses stratified Cox proportional hazards regression.
#' For reservoir subjects, it uses standard Cox regression.
#'
#' @export
InferenceSurvivalKKStratCoxIVWC = R6::R6Class("InferenceSurvivalKKStratCoxIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKStratCoxIVWC,
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

#' Stratified Cox Combined-Likelihood Compound Inference for KK Designs
#'
#' Fits the combined stratified Cox partial likelihood (matched pairs + reservoir).
#'
#' @export
InferenceSurvivalKKStratCoxOneLik = R6::R6Class("InferenceSurvivalKKStratCoxOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKStratCoxOneLik,
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
			if (should_run_asserts()) {
				assertFlag(include_covariates, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose)
			
			if (is.null(include_covariates)) {
				include_covariates = des_obj$has_covariates()
			}
			private$include_covariates_flag = include_covariates
		}
	),
	private = list(
		include_covariates_flag = NULL,
		include_covariates = function() private$include_covariates_flag
	)
)
