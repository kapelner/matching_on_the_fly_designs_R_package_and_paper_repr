#' Conditional-Logistic Inference for KK Designs with IVWC
#'
#' Fits a conditional-logistic regression for binary (incidence) responses under
#' a KK matching-on-the-fly design using the Independent-Variables-as-Working-Covariates
#' (IVWC) approach.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKClogitIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKClogitIVWC = R6::R6Class("InferenceIncidKKClogitIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitIVWC,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOneKK14} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
		}
	)
)

#' Conditional-Logistic Inference for KK Designs with Combined Likelihood
#'
#' Fits a conditional-logistic regression for binary (incidence) responses under
#' a KK matching-on-the-fly design using the combined-likelihood approach.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKClogitOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKClogitOneLik = R6::R6Class("InferenceIncidKKClogitOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitOneLik,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOneKK14} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
		}
	)
)
