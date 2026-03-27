#' Multivariate G-Computation Risk-Difference Inference for Binary Responses
#'
#' @description
#' Fits a logistic working model for an incidence outcome using treatment and all
#' recorded covariates, then estimates the marginal risk difference by
#' standardizing predicted risks under all-treated and all-control assignments
#' over the empirical covariate distribution.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignBernoulli$
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
#' infer <- DesignInferenceIncidMultiGCompRiskDiff$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferenceIncidMultiGCompRiskDiff = R6::R6Class("DesignInferenceIncidMultiGCompRiskDiff",
	inherit = DesignInferenceIncidUnivGCompRiskDiff,
	public = list(

		#' @description
		#' Initialize the multivariate g-computation RD inference object.
		#' @param des_obj A completed \code{SeqDesign} object with an incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
		}
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)

#' Multivariate G-Computation Risk-Ratio Inference for Binary Responses
#'
#' @description
#' Fits a logistic working model for an incidence outcome using treatment and all
#' recorded covariates, then estimates the marginal risk ratio by standardizing
#' predicted risks under all-treated and all-control assignments over the
#' empirical covariate distribution.
#'
#' @details
#' The point estimate is returned on the risk-ratio scale. Confidence intervals
#' and p-values use the delta method on the log-risk-ratio scale and then map
#' back to the risk-ratio scale.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignBernoulli$
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
#' infer <- DesignInferenceIncidMultiGCompRiskRatio$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferenceIncidMultiGCompRiskRatio = R6::R6Class("DesignInferenceIncidMultiGCompRiskRatio",
	inherit = DesignInferenceIncidUnivGCompRiskRatio,
	public = list(

		#' @description
		#' Initialize the multivariate g-computation RR inference object.
		#' @param des_obj A completed \code{SeqDesign} object with an incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
		}
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
