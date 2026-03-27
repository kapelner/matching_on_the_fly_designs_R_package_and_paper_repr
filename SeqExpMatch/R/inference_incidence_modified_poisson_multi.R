#' Multivariate Modified Poisson Inference for Binary Responses
#'
#' @description
#' Fits the classic modified Poisson model for binary (incidence) responses under
#' non-KK designs using treatment and all recorded covariates as predictors in a
#' Poisson log-link working model with Huber-White sandwich variance. The
#' treatment effect is reported on the log-risk-ratio scale.
#'
#' @inherit InferenceRand methods
#' @inherit InferenceBoot methods
#' @inherit InferenceAsymp methods
#' @inherit InferenceRandCI methods
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$
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
#' infer <- InferenceIncidMultiModifiedPoisson$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceIncidMultiModifiedPoisson = R6::R6Class("InferenceIncidMultiModifiedPoisson",
	inherit = InferenceIncidUnivModifiedPoisson,
	public = list(

		#' @description
		#' Initialize a multivariate modified-Poisson inference object for a
		#' completed non-KK design with a binary response.
		#' @param des_obj A completed non-KK \code{DesignSeqOneByOne} object with an
		#'   incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and
		#'   randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
		},

		#' @description
		#' Computes the modified-Poisson estimate of the treatment effect on the
		#' log-risk-ratio scale.
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
