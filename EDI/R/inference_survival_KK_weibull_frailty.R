#' Univariate Weibull Frailty IVWC Inference for KK Designs
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses using a Gamma-frail Weibull AFT model for matched pairs
#' and a standard Weibull AFT model for the reservoir. The estimates are combined via IVWC.
#' No covariates are used in the model.
#'
#' @export
InferenceSurvivalUnivKKWeibullFrailtyIVWC = R6::R6Class("InferenceSurvivalUnivKKWeibullFrailtyIVWC",
	inherit = InferenceAbstractKKWeibullFrailtyIVWC,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, verbose = FALSE){
			super$initialize(des_obj, model_formula = ~ 1, verbose = verbose)
		}
	)
)

#' Multivariate Weibull Frailty IVWC Inference for KK Designs
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses using a Gamma-frail Weibull AFT model for matched pairs
#' and a standard Weibull AFT model for the reservoir. The estimates are combined via IVWC.
#' Covariates from the design object are included in both components.
#'
#' @export
InferenceSurvivalMultiKKWeibullFrailtyIVWC = R6::R6Class("InferenceSurvivalMultiKKWeibullFrailtyIVWC",
	inherit = InferenceAbstractKKWeibullFrailtyIVWC,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	)
)

#' Univariate Weibull Frailty Combined-Likelihood Inference for KK Designs
#'
#' This class fits a single joint Weibull frailty model over all KK design data
#' (matched pairs + reservoir) for survival responses. No covariates are used.
#'
#' @export
InferenceSurvivalUnivKKWeibullFrailtyOneLik = R6::R6Class("InferenceSurvivalUnivKKWeibullFrailtyOneLik",
	inherit = InferenceAbstractKKWeibullFrailtyOneLik,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, verbose = FALSE){
			super$initialize(des_obj, model_formula = ~ 1, verbose = verbose)
		}
	)
)

#' Multivariate Weibull Frailty Combined-Likelihood Inference for KK Designs
#'
#' This class fits a single joint Weibull frailty model over all KK design data
#' (matched pairs + reservoir) for survival responses. Covariates from the design
#' object are included.
#'
#' @export
InferenceSurvivalMultiKKWeibullFrailtyOneLik = R6::R6Class("InferenceSurvivalMultiKKWeibullFrailtyOneLik",
	inherit = InferenceAbstractKKWeibullFrailtyOneLik,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	)
)
