#' GEE Inference for KK Designs with Count Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using an internal Rcpp
#' solver or \pkg{geepack}) for Poisson (count) responses under a KK 
#' matching-on-the-fly design using the treatment indicator and, optionally, 
#' all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountPoissonKKGEE$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountPoissonKKGEE = R6::R6Class("InferenceCountPoissonKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param use_rcpp Whether to use the internal Rcpp solver.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE){
			if (should_run_asserts() && !use_rcpp) {
				if (!check_package_installed("geepack")){
					stop("Package 'geepack' is required for ", class(self)[1], ". Please install it.")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, use_rcpp = use_rcpp)
		}
	),
	private = list(
		gee_response_type = function() "count",
		gee_family        = function() stats::poisson(link = "log")
	)
)
#' Univariate GEE Inference for KK Designs with Count Response
#'
#' Fits a Poisson GEE model using only the treatment indicator as a predictor.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountPoissonMultiKKGEE$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountPoissonUnivKKGEE = R6::R6Class("InferenceCountPoissonUnivKKGEE",
	lock_objects = FALSE,
	inherit = InferenceCountPoissonKKGEE,
	private = list(
		gee_predictors_df = function() data.frame(w = private$w)
	)
)
#' Multivariate GEE Inference for KK Designs with Count Response
#'
#' Fits a Poisson GEE model using the treatment indicator and all recorded 
#' covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountPoissonMultiKKGEE$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountPoissonMultiKKGEE = R6::R6Class("InferenceCountPoissonMultiKKGEE",
	lock_objects = FALSE,
	inherit = InferenceCountPoissonKKGEE
)
