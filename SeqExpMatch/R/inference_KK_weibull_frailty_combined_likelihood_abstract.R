# Abstract class for Weibull Frailty Combined-Likelihood Compound Inference (stub)
#
# @description
# Stub placeholder for a combined-likelihood compound estimator for KK matching-on-the-fly
# designs with survival responses using the Weibull AFT model.
# Methods are not yet implemented.
#
# @keywords internal
SeqDesignInferenceAbstractKKWeibullFrailtyCombinedLikelihood = R6::R6Class("SeqDesignInferenceAbstractKKWeibullFrailtyCombinedLikelihood",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		# @description
		# Initialize the inference object.
		# @param seq_des_obj		A SeqDesign object (must be a KK design).
		# @param num_cores			Number of CPU cores for parallel processing.
		# @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			if (!requireNamespace("parfm", quietly = TRUE)) {
				stop("Package 'parfm' is required for ", class(self)[1], ". Please install it.")
			}
			assertResponseType(seq_des_obj$get_response_type(), "survival")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		compute_treatment_estimate = function()
			stop(paste(class(self)[1], ": combined-likelihood estimate not yet implemented.")),
		compute_mle_confidence_interval = function(alpha = 0.05)
			stop(paste(class(self)[1], ": combined-likelihood confidence interval not yet implemented.")),
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0)
			stop(paste(class(self)[1], ": combined-likelihood p-value not yet implemented."))
	)
)
