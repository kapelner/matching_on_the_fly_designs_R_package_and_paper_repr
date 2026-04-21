#' Quantile Regression Combined-Likelihood Compound Estimator for KK Designs (Proportion)
#'
#' Fits the combined stacked quantile regression (matched-pair differences + reservoir)
#' using the treatment indicator and all recorded covariates for proportion responses.
#' Responses are transformed via logit before regression; the estimated treatment
#' effect is a log-odds-ratio shift at quantile \code{tau}.
#' Minimises the joint check-function loss over both data sources simultaneously.
#' Inference is based on the stacked combined-likelihood quantile-regression fit.
#'
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
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0.10, 0.25, 0.20, 0.40, 0.35, 0.55, 0.60, 0.75))
#' infer <- InferencePropMultiKKQuantileRegrOneLik$
#'   new(seq_des, verbose =
#' FALSE)
#' infer
#'
InferencePropMultiKKQuantileRegrOneLik = R6::R6Class("InferencePropMultiKKQuantileRegrOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKQuantileRegrOneLik,
	public = list(
		#' @description	Initialize the inference object.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param tau                             The quantile level on the logit scale, strictly
		#'   between 0 and 1. Default is 0.5.
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		initialize = function(des_obj, model_formula = NULL, tau = 0.5,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "proportion")
			}
			super$initialize(des_obj, tau, qlogis, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$y = .sanitize_proportion_response(private$y, interior = TRUE)
			if (should_run_asserts()) {
				assertNumeric(private$y, any.missing = FALSE, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			}
			# Rebuild KK summary data after sanitizing the response, otherwise the
			# superclass cache can retain raw 0/1 values and qlogis() will produce
			# non-finite values during the quantile fit and its randomization refits.
			private$cached_values$KKstats = NULL
			private$compute_basic_match_data()
		},


		#' @description
		#' Returns the estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE) super$compute_estimate()


	)
)

#' Quantile Regression Combined-Likelihood Estimator for KK Designs with Proportion Outcomes
#'
#' Public collapsed-name wrapper for \code{InferencePropMultiKKQuantileRegrOneLik}.
#'
#' @export
InferencePropKKQuantileRegrOneLik = R6::R6Class("InferencePropKKQuantileRegrOneLik",
	lock_objects = FALSE,
	inherit = InferencePropMultiKKQuantileRegrOneLik,
	public = list()
)
