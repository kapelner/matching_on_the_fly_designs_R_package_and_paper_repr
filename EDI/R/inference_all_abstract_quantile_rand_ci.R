#' Abstract mixin: Zhang combined randomisation CI for quantile regression
#'
#' Provides \code{compute_rand_confidence_interval()}
#' via Zhang's combined test-inversion method for both Bernoulli (\eqn{m = 0},
#' all subjects in the reservoir) and KK matching-on-the-fly designs
#' (\eqn{m > 0}).
#'
#' @keywords internal
InferenceAbstractQuantileRandCI = R6::R6Class("InferenceAbstractQuantileRandCI",
	lock_objects = FALSE,
	inherit = InferenceKKPassThroughCompoundNoParamBootstrap,
	public = utils::modifyList(as.list(InferenceMixinQuantileRandCI$public), list(
		#' @description Initialize the inference object.
		#' @param des_obj         A DesignSeqOneByOne object.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param verbose         Whether to print messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = NULL){
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
		},
		#' @description Computes a 1-alpha level frequentist confidence interval for the randomization test.
		#' @param alpha The significance level.
		#' @param r Number of randomization vectors.
		#' @param pval_epsilon The bisection algorithm tolerance.
		#' @param show_progress Show a text progress indicator.
		#' @param ci_search_control Optional randomization-CI search control list.
		#' @return A 1 - alpha sized frequentist confidence interval.
		compute_rand_confidence_interval = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE, ci_search_control = NULL){
			InferenceMixinQuantileRandCI$public$compute_rand_confidence_interval(alpha, r, pval_epsilon, show_progress, ci_search_control)
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinQuantileRandCI$private), list())
)
