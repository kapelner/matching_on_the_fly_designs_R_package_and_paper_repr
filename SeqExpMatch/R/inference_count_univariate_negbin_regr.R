#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceCountUnivNegBinRegr = R6::R6Class("SeqDesignInferenceCountUnivNegBinRegr",
	inherit = SeqDesignInferenceMLEorKMSummaryTable,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
        #' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
		#' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 							for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#'
		initialize = function(seq_des_obj, num_cores = 1, verbose = TRUE){			
			assertResponseType(seq_des_obj$get_response_type(), "count")			
			super$initialize(seq_des_obj, num_cores, verbose)	
			assertNoCensoring(private$any_censoring)
			private$cached_values = super$get_cached_values()
		}
	),
	
	private = list(		
		cached_values = list(),
		
		shared = function(){
			private$shared_negative_binomial_regression_inference(
				data_obj = data.frame(y = private$seq_des_obj_priv_int$y, w = private$seq_des_obj_priv_int$w)
			)
		},
		
		shared_negative_binomial_regression_inference = function(data_obj){
			tryCatch({				
				private$cached_values$summary_table = 
					coef(summary_glm_lean(suppressWarnings(MASS::glm.nb(private$seq_des_obj_priv_int$y ~ ., data = data_obj))))
			}, error = function(e){
				private$cached_values$summary_table = matrix(NA, nrow = 2, ncol = 4)				
			})
		}
		
	)		
)		