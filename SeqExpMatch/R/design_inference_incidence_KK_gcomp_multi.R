#' Multivariate G-Computation Risk-Difference Inference for KK Designs with Binary Responses
#'
#' @description
#' Fits an all-subject logistic working model for a KK incidence outcome using
#' treatment and all recorded covariates, then estimates the marginal risk
#' difference by standardizing predicted risks under all-treated and all-control
#' assignments over the empirical covariate distribution. Matched pairs are
#' treated as clusters and reservoir subjects are treated as singletons when
#' computing the sandwich covariance.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignKK14$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "incidence",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 0, 1, 0, 1, 1, 0))
#' infer <- DesignInferenceIncidMultiKKGCompRiskDiff$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferenceIncidMultiKKGCompRiskDiff = R6::R6Class("DesignInferenceIncidMultiKKGCompRiskDiff",
	inherit = DesignInferenceIncidUnivKKGCompRiskDiff,
	public = list(
		#' @description
		#' Initialize the multivariate KK g-computation RD inference object.
		#' @param seq_des_obj A completed KK \code{SeqDesign} object with an incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)

#' Multivariate G-Computation Risk-Ratio Inference for KK Designs with Binary Responses
#'
#' @description
#' Fits an all-subject logistic working model for a KK incidence outcome using
#' treatment and all recorded covariates, then estimates the marginal risk ratio
#' by standardizing predicted risks under all-treated and all-control assignments
#' over the empirical covariate distribution. Matched pairs are treated as
#' clusters and reservoir subjects are treated as singletons when computing the
#' sandwich covariance.
#'
#' @details
#' The point estimate is returned on the risk-ratio scale. Confidence intervals
#' and p-values use the delta method on the log-risk-ratio scale and then map
#' back to the risk-ratio scale.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignKK14$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "incidence",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 0, 1, 0, 1, 1, 0))
#' infer <- DesignInferenceIncidMultiKKGCompRiskRatio$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferenceIncidMultiKKGCompRiskRatio = R6::R6Class("DesignInferenceIncidMultiKKGCompRiskRatio",
	inherit = DesignInferenceIncidUnivKKGCompRiskRatio,
	public = list(
		#' @description
		#' Initialize the multivariate KK g-computation RR inference object.
		#' @param seq_des_obj A completed KK \code{SeqDesign} object with an incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
