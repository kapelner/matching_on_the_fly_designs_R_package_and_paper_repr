#' G-Computation Risk-Difference Inference for KK Designs with Binary Responses
#'
#' Fits an all-subject logistic working model for a KK incidence outcome using
#' treatment and, optionally, all recorded covariates, then estimates the marginal
#' risk difference by standardizing predicted risks under all-treated and
#' all-control assignments over the empirical covariate distribution. Matched
#' pairs are treated as clusters and reservoir subjects are treated as singletons
#' when computing the sandwich covariance.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKGCompRiskDiff$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKGCompRiskDiff = R6::R6Class("InferenceIncidKKGCompRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceIncidKKGCompAbstract,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		},

		get_estimand_type = function() "RD"
	)
)

#' G-Computation Risk-Ratio Inference for KK Designs with Binary Responses
#'
#' Fits a all-subject logistic working model for a KK incidence outcome using
#' treatment and, optionally, all recorded covariates, then estimates the marginal
#' risk ratio by standardizing predicted risks under all-treated and all-control
#' assignments over the empirical covariate distribution. Matched pairs are
#' treated as clusters and reservoir subjects are treated as singletons when
#' computing the sandwich covariance.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKGCompRiskRatio$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKGCompRiskRatio = R6::R6Class("InferenceIncidKKGCompRiskRatio",
	lock_objects = FALSE,
	inherit = InferenceIncidKKGCompAbstract,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		},

		get_estimand_type = function() "RR"
	)
)

#' Modified-Poisson Inference for KK Designs with Binary Responses
#'
#' Fits an all-subject modified-Poisson working model for incidence outcomes under
#' a KK matching-on-the-fly design using treatment and, optionally, all recorded
#' covariates as predictors. Matched pairs are treated as clusters and reservoir
#' subjects are treated as singleton clusters when computing the sandwich
#' covariance.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKModifiedPoisson$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKModifiedPoisson = R6::R6Class("InferenceIncidKKModifiedPoisson",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKModifiedPoisson,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
