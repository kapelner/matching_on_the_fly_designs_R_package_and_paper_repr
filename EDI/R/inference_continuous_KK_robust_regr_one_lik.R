#' Robust-Regression Combined-Likelihood Inference for KK Designs
#'
#' Fits a single stacked robust regression over matched-pair differences and
#' reservoir observations for KK matching-on-the-fly designs with continuous
#' responses, using the treatment indicator and, optionally, all recorded
#' covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceContinKKRobustRegrOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceContinKKRobustRegrOneLik = R6::R6Class("InferenceContinKKRobustRegrOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKRobustRegrOneLik,
	public = list(
	)
)
