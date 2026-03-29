#' Multivariate Modified-Poisson Inference for KK Designs with Binary Responses
#'
#' Fits an all-subject modified-Poisson working model for incidence outcomes under
#' a KK matching-on-the-fly design using treatment and all recorded covariates as
#' predictors. Matched pairs are treated as clusters and reservoir subjects are
#' treated as singleton clusters when computing the sandwich covariance, so the
#' estimated treatment effect is a log risk ratio with cluster-robust inference.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneKK14$
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
#' infer <- InferenceIncidMultiKKModifiedPoisson$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceIncidMultiKKModifiedPoisson = R6::R6Class("InferenceIncidMultiKKModifiedPoisson",
	lock_objects = FALSE,
	inherit = InferenceIncidUnivKKModifiedPoisson,
	public = list(

	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
