# Abstract Quantile Regression Combined-Likelihood Compound Estimator for KK Designs (stub)
#
# @description
# Stub placeholder for a combined-likelihood compound quantile regression estimator for KK
# matching-on-the-fly designs. Methods are not yet implemented.
#
# @keywords internal
SeqDesignInferenceAbstractKKQuantileRegrCombinedLikelihood = R6::R6Class("SeqDesignInferenceAbstractKKQuantileRegrCombinedLikelihood",
	inherit = SeqDesignInferenceAbstractQuantileRandCI,
	public = list(

		# @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		# @param tau				The quantile level for regression, strictly between 0 and 1.
		# @param transform_y_fn	A function applied to y values before quantile regression.
		# @param num_cores			The number of CPU cores to use to parallelize sampling.
		# @param verbose			A flag indicating whether messages should be displayed. Default is FALSE.
		initialize = function(seq_des_obj, tau = 0.5, transform_y_fn = identity, num_cores = 1, verbose = FALSE){
			assertNumeric(tau, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			if (!requireNamespace("quantreg", quietly = TRUE)) {
				stop("Package 'quantreg' is required. Please install it with install.packages(\"quantreg\").")
			}
			private$tau = tau
			private$transform_y_fn_list = list(fn = transform_y_fn)
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		compute_treatment_estimate = function()
			stop(paste(class(self)[1], ": combined-likelihood estimate not yet implemented.")),
		compute_mle_confidence_interval = function(alpha = 0.05)
			stop(paste(class(self)[1], ": combined-likelihood confidence interval not yet implemented.")),
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0)
			stop(paste(class(self)[1], ": combined-likelihood p-value not yet implemented."))
	),

	private = list(
		tau = NULL,
		transform_y_fn_list = NULL
	)
)
