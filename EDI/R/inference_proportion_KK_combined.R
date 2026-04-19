#' GEE Inference for KK Designs with Proportion Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for proportion (continuous values in (0, 1)) responses under a KK
#' matching-on-the-fly design using the treatment indicator and, optionally,
#' all recorded covariates as predictors.
#'
#' @details
#' This class requires the \pkg{geepack} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{geepack} before using this class.
#'
#' @export
InferencePropKKGEE = R6::R6Class("InferencePropKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				if (!check_package_installed("geepack")){
					stop("Package 'geepack' is required for ", class(self)[1], ". Please install it.")
				}
			}
			super$initialize(des_obj, include_covariates, verbose)
		}
	),
	private = list(
		gee_response_type = function() "proportion",
		gee_family        = function() stats::binomial(link = "logit")
	)
)

#' GLMM Inference for KK Designs with Proportion Response
#'
#' Fits a Generalized Linear Mixed Model (GLMM) using the \pkg{glmmTMB} fitter for
#' proportion (continuous values in (0, 1)) responses under a KK
#' matching-on-the-fly design using the treatment indicator and, optionally, all
#' recorded covariates as fixed-effect predictors.
#'
#' @details
#' This class requires the \pkg{glmmTMB} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{glmmTMB} before using this class.
#'
#' @export
InferencePropKKGLMM = R6::R6Class("InferencePropKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			super$initialize(des_obj, include_covariates, verbose)
		}
	),
	private = list(
		glmm_response_type  = function() "proportion",
		glmm_family         = function() stats::binomial(link = "logit")
	)
)
