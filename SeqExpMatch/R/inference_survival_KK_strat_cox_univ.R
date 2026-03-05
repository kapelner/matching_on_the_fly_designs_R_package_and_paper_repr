#' Univariate Stratified Cox / Standard Cox Compound Inference for KK Designs
#'
#' @description
#' Fits a compound estimator for KK matching-on-the-fly designs with survival responses
#' using only the treatment indicator. For matched pairs, it uses stratified Cox
#' proportional hazards regression (each pair is a stratum). For reservoir subjects, it
#' uses standard Cox regression. Both models estimate log-hazard ratios, which are then
#' combined via a variance-weighted linear combination.
#'
#' @export
SeqDesignInferenceSurvivalUnivKKStratCox = R6::R6Class("SeqDesignInferenceSurvivalUnivKKStratCox",
	inherit = SeqDesignInferenceAbstractKKStratCox,
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
