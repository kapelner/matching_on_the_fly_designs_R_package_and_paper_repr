#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceCountUnivNegBinRegr = R6::R6Class("SeqDesignInferenceCountUnivNegBinRegr",
	inherit = SeqDesignInferenceMLEorKMforGLMs,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference 
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs), 
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "count")			
			super$initialize(seq_des_obj, num_cores, verbose)
		},
		
		#' @description
		#' Computes the appropriate estimate
		#' 
		#' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
		compute_treatment_estimate = function(){
			private$shared() # Ensure the model is fitted and cached values are populated
			private$cached_values$beta_hat_T
		}
	),
	
	private = list(		
		generate_mod = function(){
			fast_negbin_regression_with_var(
				Xmm = cbind(1, private$w), # just Intercept + Treatment
				y = private$y
			)
		}		
	)		
)