#' Internal Base Class for KK Matching-on-the-Fly Designs
#'
#' @name InferenceKKPassThroughCompound
#' @description Internal method.
#' An abstract R6 class that provides relevant methods when the designs are KK matching-on-the-fly.
#'
#' @keywords internal
InferenceKKPassThroughCompound = R6::R6Class("InferenceKKPassThroughCompound",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = c(InferenceMixinKKPassThrough$public, InferenceMixinKKPassThroughCompound$public, list(
		#' @description Initialize
		#' @param des_obj         A DesignSeqOneByOne object.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param verbose         Whether to print messages.
		#' @param harden          Whether to apply robustness measures.
		#' @param smart_default   Whether to use smart optimizer start values.
		initialize = function(des_obj, verbose = FALSE, harden = TRUE, model_formula = NULL, smart_default = TRUE){
			super$initialize(des_obj, verbose = verbose, harden = harden, model_formula = model_formula, smart_default = smart_default)
			private$init_kk_passthrough(des_obj)
		}
	)),
	private = c(InferenceMixinKKPassThrough$private, InferenceMixinKKPassThroughCompound$private, list())
)
