#' Partial Proportional Odds Regression Inference for Ordinal Responses
#'
#' Fits a partial proportional odds regression for ordinal responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalPartialProportionalOddsRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalPartialProportionalOddsRegr = R6::R6Class("InferenceOrdinalPartialProportionalOddsRegr",
	lock_objects = FALSE,
	inherit = InferenceOrdinalPartialProportionalOddsAbstract,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param nonparallel Covariate names that may vary across thresholds.
		#' @param verbose Whether to print progress messages.
		#' @param harden Whether to apply robustness measures.
		#' @param smart_default Whether to use smart optimizer starts.
		initialize = function(des_obj, verbose = FALSE, harden = TRUE, model_formula = NULL, nonparallel = character(0), smart_default = TRUE){
			super$initialize(des_obj, verbose = verbose, harden = harden, model_formula = model_formula, nonparallel = nonparallel, smart_default = smart_default)
		}
	),

	private = list(
		ppo_covariate_matrix = function(){
			private$X
		},

		build_design_matrix = function(){
			X_cov = private$ppo_covariate_matrix()
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = matrix(private$w, ncol = 1L)
				colnames(X) = "treatment"
			} else {
				X = cbind(treatment = private$w, X_cov)
			}
			X
		}
	)
)
