#' LWA-style Marginal Cox IVWC Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with survival responses
#' using a marginal Cox model with Lee-Wei-Amato style cluster-robust variance for
#' matched pairs and standard Cox regression for reservoir subjects.
#'
#' @export
InferenceSurvivalKKLWACoxIVWC = R6::R6Class("InferenceSurvivalKKLWACoxIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKLWACoxIVWC,
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

#' LWA-style Marginal Cox Combined-Likelihood Inference for KK Designs
#'
#' Fits a combined-likelihood Cox model for KK matching-on-the-fly designs with
#' survival responses using a marginal approach over all subjects.
#'
#' @export
InferenceSurvivalKKLWACoxOneLik = R6::R6Class("InferenceSurvivalKKLWACoxOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKLWACoxOneLik,
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
