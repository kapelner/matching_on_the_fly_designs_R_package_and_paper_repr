#' Multivariate Survival Rank-based Regression (AFT) Compound Inference for KK Designs
#'
#' @description
#' Fits a robust compound estimator for KK matching-on-the-fly designs with survival
#' responses using rank-based estimating equations with the treatment indicator and
#' all recorded covariates. For matched pairs, it uses a rank-based AFT model with
#' clustering. For reservoir subjects, it uses a standard rank-based AFT model. Both
#' models estimate log-time ratios, which are then combined via a variance-weighted
#' linear combination. This is the multivariate (covariate-adjusted) variant.
#'
#' @export
SeqDesignInferenceSurvivalMultiKKRankRegrIVWC = R6::R6Class("SeqDesignInferenceSurvivalMultiKKRankRegrIVWC",
	inherit = SeqDesignInferenceAbstractKKSurvivalRankRegrIVWC,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design).
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),
	private = list(
		include_covariates = function() TRUE
	)
)
