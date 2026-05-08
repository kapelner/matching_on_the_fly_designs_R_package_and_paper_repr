#' Robust-Regression IVWC Compound Inference for KK Designs
#'
#' Fits a variance-weighted compound estimator for KK matching-on-the-fly designs
#' with continuous responses using robust regression for matched-pair differences
#' and reservoir outcomes, with treatment and, optionally, all recorded covariates
#' as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceContinKKRobustRegrIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceContinKKRobustRegrIVWC = R6::R6Class("InferenceContinKKRobustRegrIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKRobustRegrIVWC,
	public = list(
	)
)
