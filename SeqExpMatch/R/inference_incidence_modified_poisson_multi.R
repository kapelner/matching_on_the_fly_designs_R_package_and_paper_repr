#' Multivariate Modified Poisson Inference for Binary Responses
#'
#' Fits the classic modified Poisson model for binary (incidence) responses under
#' non-KK designs using treatment and all recorded covariates as predictors in a
#' Poisson log-link working model with Huber-White sandwich variance. The
#' treatment effect is reported on the log-risk-ratio scale.
#'
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
	lock_objects = FALSE,
	inherit = InferenceIncidUnivModifiedPoisson,
	public = list(


	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
