#' Univariate Ordered Probit Inference for Ordinal Responses
#'
#' @description
#' Ordinal probit (ordered probit) model inference for ordinal responses.
#'
#' @export
SeqDesignInferenceOrdinalUniOrderedProbitRegr = R6::R6Class("SeqDesignInferenceOrdinalUniOrderedProbitRegr",
	inherit = SeqDesignInferenceOrdinalUniCumulProbitRegr,
	public = list(
		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param seq_des_obj A SeqDesign object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param num_cores The number of CPU cores to use.
		#' @param verbose A flag indicating whether messages should be displayed.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	)
)

#' Multivariate Ordered Probit Inference for Ordinal Responses
#'
#' @description
#' Ordinal probit (ordered probit) inference for ordinal responses with
#' treatment and observed covariates entering linearly into the latent normal
#' index.
#'
#' @export
SeqDesignInferenceOrdinalMultiOrderedProbitRegr = R6::R6Class("SeqDesignInferenceOrdinalMultiOrderedProbitRegr",
	inherit = SeqDesignInferenceOrdinalMultiCumulProbitRegr,
	public = list(
		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param seq_des_obj A SeqDesign object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param num_cores The number of CPU cores to use.
		#' @param verbose A flag indicating whether messages should be displayed.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	)
)
