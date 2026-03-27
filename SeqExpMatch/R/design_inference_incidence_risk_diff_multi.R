#' Multivariate Risk-Difference Regression for Binary Responses
#'
#' @description
#' Fits a direct risk-difference estimator for binary (incidence) responses
#' under non-KK designs using treatment and all recorded covariates in a linear
#' probability model with HC2 heteroskedasticity-robust variance. The treatment
#' effect is reported on the risk-difference scale.
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
#' infer <- DesignInferenceIncidMultiRiskDiff$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferenceIncidMultiRiskDiff = R6::R6Class("DesignInferenceIncidMultiRiskDiff",
	inherit = DesignInferenceIncidUnivRiskDiff,
	public = list(

		#' @description
		#' Initialize a multivariate risk-difference inference object for a
		#' completed non-KK design with a binary response.
		#' @param des_obj A completed non-KK \code{SeqDesign} object with an
		#'   incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and
		#'   randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
		},

		#' @description
		#' Computes the direct risk-difference estimate of the treatment effect.
		compute_treatment_estimate = function(){
			super$compute_treatment_estimate()
		}
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
