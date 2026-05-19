#' Internal Base Class for KK Matching-on-the-Fly Designs
#'
#' @name InferenceKKPassThroughCompound
#' @description Internal method.
#' An abstract R6 class that provides relevant methods when the designs are KK matching-on-the-fly.
#'
#' @keywords internal
InferenceKKPassThroughCompound = R6::R6Class("InferenceKKPassThroughCompound",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = utils::modifyList(utils::modifyList(as.list(InferenceMixinKKPassThrough$public), as.list(InferenceMixinKKPassThroughCompound$public)), list(
		#' @description Initialize
		#' @param des_obj         A DesignSeqOneByOne object.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param verbose         Whether to print messages.
		#' @param harden          Whether to apply robustness measures.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, verbose = FALSE, harden = TRUE, model_formula = NULL, smart_cold_start_default = TRUE){
			super$initialize(des_obj, verbose = verbose, harden = harden, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKPassThrough$private), as.list(InferenceMixinKKPassThroughCompound$private))
)
