#' Conditional Logistic Regression IVWC Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with binary (incidence)
#' responses using the treatment indicator and, optionally, all recorded covariates.
#' For matched pairs, a conditional logistic regression model is used (via the internal
#' \code{clogit_helper}). For reservoir subjects, a standard logistic regression
#' is used. The two estimates are combined via a variance-weighted linear combination.
#'
#' @export
InferenceIncidKKClogitIVWC = R6::R6Class("InferenceIncidKKClogitIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitIVWC,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
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

#' Conditional Logistic Combined-Likelihood Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with binary (incidence)
#' responses using the treatment indicator and, optionally, all recorded covariates.
#' Uses the combined logistic likelihood over discordant matched-pair differences
#' and reservoir subjects.
#'
#' @export
InferenceIncidKKClogitOneLik = R6::R6Class("InferenceIncidKKClogitOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitOneLik,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
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
