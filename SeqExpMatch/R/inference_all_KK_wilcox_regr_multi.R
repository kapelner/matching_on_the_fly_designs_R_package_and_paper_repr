#' Multivariate Wilcox Rank-based Regression Compound Inference for KK Designs
#'
#' @description
#' Fits a robust compound estimator for KK matching-on-the-fly designs using rank-based
#' regression (R-estimation) with the treatment indicator and all recorded covariates.
#' For matched pairs, it uses R-estimation on the within-pair differences of responses
#' and covariates. For reservoir subjects, it uses a standard rank-based linear model
#' including both treatment and covariates. This is the multivariate (covariate-adjusted)
#' variant.
#'
#' @export
SeqDesignInferenceAllKKWilcoxRegrMulti = R6::R6Class("SeqDesignInferenceAllKKWilcoxRegrMulti",
	inherit = SeqDesignInferenceAbstractKKWilcoxRegr,
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
		include_covariates = function() TRUE
	)
)
