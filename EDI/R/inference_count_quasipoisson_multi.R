#' Multivariate Quasi-Poisson Regression Inference for Count Responses
#'
#' Fits a quasi-Poisson log-link regression for count responses using the treatment
#' indicator and all recorded covariates. The treatment effect is reported on the
#' log-rate scale and inference uses the model-based quasi-Poisson variance with an
#' estimated dispersion parameter.
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
#'   response_type = "count",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 1, 2, 2, 3, 3, 4))
#' infer <- InferenceCountMultiQuasiPoissonRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceCountMultiQuasiPoissonRegr = R6::R6Class("InferenceCountMultiQuasiPoissonRegr",
	lock_objects = FALSE,
	inherit = InferenceCountUnivQuasiPoissonRegr,
	public = list(

	),

	private = list(
		build_design_matrix = function(){
			Xmm = private$create_design_matrix()
			full_names = c("(Intercept)", "treatment", if (ncol(Xmm) > 2L) paste0("x", seq_len(ncol(Xmm) - 2L)) else NULL)
			colnames(Xmm) = full_names[seq_len(ncol(Xmm))]
			Xmm
		}
	)
)
