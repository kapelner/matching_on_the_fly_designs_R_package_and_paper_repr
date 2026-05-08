#' Stratified Cox / Standard Cox Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with survival responses.
#' For matched pairs, it uses stratified Cox proportional hazards regression.
#' For reservoir subjects, it uses standard Cox regression.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalKKStratCoxIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalKKStratCoxIVWC = R6::R6Class("InferenceSurvivalKKStratCoxIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKStratCoxIVWC,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	)
)

#' Stratified Cox Combined-Likelihood Compound Inference for KK Designs
#'
#' Fits the combined stratified Cox partial likelihood (matched pairs + reservoir).
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalKKStratCoxOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalKKStratCoxOneLik = R6::R6Class("InferenceSurvivalKKStratCoxOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKStratCoxOneLik,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
			}
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	)
)
