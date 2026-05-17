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
	inherit = InferenceAsymp,
	public = utils::modifyList(as.list(InferenceMixinKKGEEShared$public), list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param use_rcpp Whether to use the internal Rcpp solver.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, smart_cold_start_default = TRUE){
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_gee_shared(des_obj, use_rcpp = use_rcpp, model_formula = model_formula)
		},
		#' @description Compute the treatment estimate.
		#' @param estimate_only Whether to skip standard-error calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_gee_dispatch(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared_gee_dispatch(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Row weights for the bootstrap sample.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			private$shared_combined_bootstrap(subject_or_block_weights, estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			super$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug, bootstrap_type)
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKGEEShared$private), list(
		gee_response_type = function() "count",
		gee_family        = function() stats::poisson(link = "log"),
		shared_gee_dispatch = function(estimate_only = FALSE) private$shared_gee_default(estimate_only)
	))
)
#' GEE Inference for KK Designs with Multi-Count Response
#' @export
InferenceCountPoissonMultiKKGEE = R6::R6Class("InferenceCountPoissonMultiKKGEE",
	lock_objects = FALSE,
	inherit = InferenceCountPoissonKKGEE
)
#' GEE Inference for KK Designs with Univariate Count Response
#' @export
InferenceCountPoissonUnivKKGEE = R6::R6Class("InferenceCountPoissonUnivKKGEE",
	lock_objects = FALSE,
	inherit = InferenceCountPoissonKKGEE,
	private = list(
		shared_gee_dispatch = function(estimate_only = FALSE) private$shared_gee_univariate(estimate_only)
	)
)
