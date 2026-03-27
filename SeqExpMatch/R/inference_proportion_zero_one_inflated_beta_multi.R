#' Multivariate Zero/One-Inflated Beta Inference for Proportion Responses
#'
#' @description
#' Fits a zero/one-inflated beta regression for proportion responses using the
#' treatment indicator and all recorded covariates in the beta mean submodel.
#' Exact 0 and exact 1 values are modeled with separate inflation masses, and
#' interior values are modeled with a beta distribution. The reported treatment
#' effect is the treatment coefficient from the beta mean submodel, on the logit
#' scale.
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
#'   response_type = "proportion",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0.10, 0.25, 0.20, 0.40, 0.35, 0.55, 0.60, 0.75))
#' infer <- InferencePropMultiZeroOneInflatedBetaRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferencePropMultiZeroOneInflatedBetaRegr = R6::R6Class("InferencePropMultiZeroOneInflatedBetaRegr",
	inherit = InferencePropZeroOneInflatedBetaAbstract,
	public = list(

	),

	private = list(
		predictors_df = function(){
			private$get_X()
		}
	)
)
