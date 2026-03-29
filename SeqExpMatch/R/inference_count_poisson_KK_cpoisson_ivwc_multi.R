#' Multivariate Conditional Poisson / Negative Binomial Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with count responses
#' using the treatment indicator and all recorded covariates. For matched pairs, it uses
#' conditional Poisson regression. For reservoir subjects, it uses Negative Binomial regression.
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
#'   response_type = "count",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 1, 2, 2, 3, 3, 4))
#' infer <- InferenceCountPoissonMultiKKCPoissonIVWC$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceCountPoissonMultiKKCPoissonIVWC = R6::R6Class("InferenceCountPoissonMultiKKCPoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKPoissonCPoissonIVWC,
	public = list(




	),
	private = list(
		include_covariates = function() TRUE
	)
)
