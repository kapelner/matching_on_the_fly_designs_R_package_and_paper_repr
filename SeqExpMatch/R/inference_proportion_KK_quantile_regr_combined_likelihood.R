#' Quantile Regression Combined-Likelihood Compound Estimator for KK Designs (Proportion)
#'
#' @description
#' Fits the combined stacked quantile regression (matched-pair differences + reservoir)
#' using the treatment indicator and all recorded covariates for proportion responses.
#' Responses are transformed via logit before regression; the estimated treatment
#' effect is a log-odds-ratio shift at quantile \code{tau}.
#' Minimises the joint check-function loss over both data sources simultaneously.
#' Inference is based on the stacked combined-likelihood quantile-regression fit.
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
#'   response_type = "proportion",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0.10, 0.25, 0.20, 0.40, 0.35, 0.55, 0.60, 0.75))
#' infer <- InferencePropMultiKKQuantileRegrCombinedLikelihood$
#'   new(seq_des, verbose =
#' FALSE)
#' infer
#'
InferencePropMultiKKQuantileRegrCombinedLikelihood = R6::R6Class("InferencePropMultiKKQuantileRegrCombinedLikelihood",
	inherit = InferenceAbstractKKQuantileRegrCombinedLikelihood,
	public = list(
		#' @description	Initialize the inference object.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param tau                             The quantile level on the logit scale, strictly
		#'   between 0 and 1. Default is 0.5.
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(des_obj, tau = 0.5, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "proportion")
			super$initialize(des_obj, tau, qlogis, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
			assertNumeric(private$y, any.missing = FALSE, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
		},


		#' @description
		#' Returns the estimated treatment effect.
		compute_treatment_estimate = function(estimate_only = FALSE) super$compute_treatment_estimate()


	)
)
