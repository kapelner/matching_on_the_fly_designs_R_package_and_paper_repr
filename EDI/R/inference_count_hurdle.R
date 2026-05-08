#' Hurdle Poisson Regression Inference for Count Responses
#'
#' Fits a hurdle Poisson regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountHurdleNegBin$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountHurdlePoisson = R6::R6Class("InferenceCountHurdlePoisson",
	lock_objects = FALSE,
	inherit = InferenceCountZeroAugmentedPoissonAbstract,
	public = list(
	),
	private = list(
		za_family = function() glmmTMB::truncated_poisson(link = "log"),
		za_description = function() "Hurdle Poisson"
	)
)

#' Hurdle Negative Binomial Regression Inference for Count Responses
#'
#' Fits a hurdle negative binomial regression for count responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountHurdleNegBin$new(seq_des, model_formula = ~ x1)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountHurdleNegBin = R6::R6Class("InferenceCountHurdleNegBin",
	lock_objects = FALSE,
	inherit = InferenceCountHurdleNegBinAbstract,
	public = list(
	),
	private = list(
		hurdle_description = function() "Hurdle Negative Binomial"
	)
)
