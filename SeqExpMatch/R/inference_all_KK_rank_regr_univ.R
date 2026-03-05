#' Univariate Rank-based Regression (R-Estimation) Compound Inference for KK Designs
#'
#' @description
#' Fits a robust compound estimator for KK matching-on-the-fly designs using rank-based
#' regression (R-estimation) with only the treatment indicator. For matched pairs, it
#' uses R-estimation on the within-pair response differences. For reservoir subjects,
#' it uses a standard rank-based linear model with the treatment assignment as the
#' sole predictor.
#'
#' @export
SeqDesignInferenceAllKKRankRegrUniv = R6::R6Class("SeqDesignInferenceAllKKRankRegrUniv",
	inherit = SeqDesignInferenceAbstractKKRankRegr,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param seq_des_obj		A SeqDesign object (must be a KK design).
		#' @param num_cores			Number of CPU cores for parallel processing.
		#' @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),
	private = list(
		include_covariates = function() FALSE
	)
)
