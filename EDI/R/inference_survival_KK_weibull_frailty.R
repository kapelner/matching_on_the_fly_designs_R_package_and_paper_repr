#' Weibull Frailty IVWC Inference for KK Designs
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses using a Gamma-frail Weibull AFT model for matched pairs
#' and a standard Weibull AFT model for the reservoir. The estimates are combined via IVWC.
#'
#' @export
InferenceSurvivalKKWeibullFrailtyIVWC = R6::R6Class("InferenceSurvivalKKWeibullFrailtyIVWC",
	inherit = InferenceAbstractKKWeibullFrailtyIVWC,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param model_formula Optional formula for covariate adjustment. If \code{NULL} (default), 
		#'   covariates from the design object are included. Use \code{~ 1} for univariate.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	)
)

#' Weibull Frailty Combined-Likelihood Inference for KK Designs
#'
#' This class fits a single joint Weibull frailty model over all KK design data
#' (matched pairs + reservoir) for survival responses.
#'
#' @export
InferenceSurvivalKKWeibullFrailtyOneLik = R6::R6Class("InferenceSurvivalKKWeibullFrailtyOneLik",
	inherit = InferenceAbstractKKWeibullFrailtyOneLik,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param model_formula Optional formula for covariate adjustment. If \code{NULL} (default), 
		#'   covariates from the design object are included. Use \code{~ 1} for univariate.
		#' @param use_rcpp Whether to use the custom Rcpp likelihood optimizer.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg The optimization algorithm to use.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, optimization_alg = "lbfgs"){
			super$initialize(des_obj, model_formula = model_formula, use_rcpp = use_rcpp, verbose = verbose, optimization_alg = optimization_alg)
		}
	)
)
