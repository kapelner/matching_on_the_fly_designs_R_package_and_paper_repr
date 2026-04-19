#' GEE Inference for KK Designs with Binary Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for binary (incidence) responses under a KK matching-on-the-fly design using
#' the treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceIncidKKGEE = R6::R6Class("InferenceIncidKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
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
				if (!check_package_installed("geepack")){
					stop("Package 'geepack' is required for ", class(self)[1], ". Please install it.")
				}
			}
			super$initialize(des_obj, include_covariates, verbose)
		}
	),
	private = list(
		gee_response_type = function() "incidence",
		gee_family        = function() stats::binomial(link = "logit")
	)
)

#' GLMM Inference for KK Designs with Binary Response
#'
#' Fits a Generalized Linear Mixed Model (GLMM) using the \pkg{glmmTMB} fitter for
#' binary (incidence) responses under a KK matching-on-the-fly design using the
#' treatment indicator and, optionally, all recorded covariates as fixed-effect
#' predictors.
#'
#' @export
InferenceIncidKKGLMM = R6::R6Class("InferenceIncidKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
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
			super$initialize(des_obj, include_covariates, verbose)
		}
	),
	private = list(
		glmm_response_type  = function() "incidence",
		glmm_family         = function() stats::binomial(link = "logit")
	)
)
