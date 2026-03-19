#' Univariate Ordered Probit Inference for Ordinal Responses
#'
#' @description
#' Ordinal probit (ordered probit) model inference for ordinal responses.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignCRD$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- SeqDesignInferenceOrdinalUniOrderedProbitRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
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
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignCRD$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- SeqDesignInferenceOrdinalMultiOrderedProbitRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
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
