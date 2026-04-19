#' KK Hurdle Poisson IVWC Inference for Count Responses
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with count
#' responses using a hurdle-Poisson mixed model for matched pairs and an ordinary
#' Poisson regression for reservoir subjects. Estimates are combined via
#' inverse-variance weighting.
#'
#' @export
InferenceCountKKHurdlePoissonIVWC = R6::R6Class("InferenceCountKKHurdlePoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKHurdlePoissonIVWC,
	public = list(
	)
)

#' KK Poisson-Conditional-Poisson IVWC Inference for Count Responses
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with count
#' responses using a conditional Poisson model for matched pairs and an ordinary
#' Poisson regression for reservoir subjects. Estimates are combined via
#' inverse-variance weighting.
#'
#' @export
InferenceCountKKCPoissonIVWC = R6::R6Class("InferenceCountKKCPoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKPoissonCPoissonIVWC,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
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
