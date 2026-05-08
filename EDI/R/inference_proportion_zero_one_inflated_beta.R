#' Zero/One-Inflated Beta Inference for Proportion Responses
#'
#' Fits a zero/one-inflated beta regression for proportion responses using the
#' treatment indicator and, optionally, all recorded covariates in the beta mean
#' submodel. Exact 0 and exact 1 values are modeled with separate inflation
#' masses, and interior values are modeled with a beta distribution. The reported
#' treatment effect is the treatment coefficient from the beta mean submodel, on
#' the logit scale.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'proportion')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferencePropZeroOneInflatedBetaRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferencePropZeroOneInflatedBetaRegr = R6::R6Class("InferencePropZeroOneInflatedBetaRegr",
	lock_objects = FALSE,
	inherit = InferencePropZeroOneInflatedBetaAbstract,
	public = list(
	)
)
