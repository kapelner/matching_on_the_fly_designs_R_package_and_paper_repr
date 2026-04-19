#' Linear Mixed Model Inference for KK Designs with Continuous Response
#'
#' Fits a linear mixed model using the \pkg{glmmTMB} fitter for continuous responses
#' under a KK matching-on-the-fly design using the treatment indicator and, optionally,
#' all recorded covariates as fixed-effect predictors. A Gaussian identity-link working
#' model is used. The matched-pair strata enter the model as a subject-level random
#' intercept \code{(1 | group_id)}, which accounts for within-pair correlation.
#'
#' @details
#' This class requires the \pkg{glmmTMB} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{glmmTMB} before using this class.
#'
#' @export
InferenceContinKKGLMM = R6::R6Class("InferenceContinKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(

		#' @description
		#' Initialize a KK GLMM inference object.
		#' @param des_obj A completed \code{Design} object with a continuous response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as fixed-effect predictors. If \code{FALSE}, only the treatment
		#'   indicator is used. If \code{NULL} (default), it is set to \code{TRUE} if the
		#'   design contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertFlag(include_covariates, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose)
			
			if (is.null(include_covariates)) {
				include_covariates = des_obj$has_covariates()
			}
			private$include_covariates = include_covariates
		}
	),
	private = list(
		include_covariates = NULL,
		glmm_response_type = function() "continuous",
		glmm_family        = function() stats::gaussian(link = "identity"),
		
		glmm_predictors_df = function(){
			if (private$include_covariates) {
				# Use default multivariate logic from InferenceAbstractKKGLMM
				super$glmm_predictors_df()
			} else {
				# Univariate logic
				data.frame(w = private$w)
			}
		}
	)
)
