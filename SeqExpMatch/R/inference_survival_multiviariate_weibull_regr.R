#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceSurvivalMultiWeibullRegr = R6::R6Class("SeqDesignInferenceSurvivalMultiWeibullRegr",
	inherit = SeqDesignInferenceSurvivalUniWeibullRegr,
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
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),
	
	private = list(		
		shared = function(){
			private$compute_survival_weibull_regression(
				data_obj = cbind(private$seq_des_obj_priv_int$w, private$get_X())
			)
		}		
	)		
)
		