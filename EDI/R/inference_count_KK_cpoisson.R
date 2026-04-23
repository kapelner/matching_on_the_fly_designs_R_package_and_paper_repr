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
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOneKK14} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp
		#'   implementations where available. If \code{FALSE}, use \pkg{glmmTMB} for
		#'   the matched-pair component.
		#' @param optimization_alg Optimization algorithm. Default is dispatched via policy.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = NULL, verbose = FALSE){
			super$initialize(des_obj, model_formula = model_formula, use_rcpp = use_rcpp, optimization_alg = optimization_alg, verbose = verbose)
		}
	)
)

#' Conditional-Poisson Inference for KK Designs with IVWC
#'
#' Fits a conditional-Poisson regression for count responses under a KK design
#' using the Independent-Variables-as-Working-Covariates (IVWC) approach.
#'
#' @export
InferenceCountKKCPoissonIVWC = R6::R6Class("InferenceCountKKCPoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKPoissonCPoissonIVWC,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOneKK14} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	)
)

#' Conditional-Poisson Inference for KK Designs with Combined Likelihood
#'
#' Fits a conditional-Poisson regression for count responses under a KK design
#' using the combined-likelihood approach.
#'
#' @export
InferenceCountKKCPoissonOneLik = R6::R6Class("InferenceCountKKCPoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKPoissonCPoissonOneLik,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOneKK14} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	)
)
