#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceIncidUnivLogRegr = R6::R6Class("SeqDesignInferenceIncidUnivLogRegr",
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
			assertResponseType(seq_des_obj$get_response_type(), "incidence")			
			super$initialize(seq_des_obj, num_cores, verbose)	
			assertNoCensoring(private$any_censoring)
			private$cached_values = super$get_cached_values() #R6 is dumb
		}
	),
	
	private = list(		
		cached_values = list(),
		
		shared = function(){
			private$cached_values$summary_table = tryCatch({
					private$compute_summary_table()
				}, error = function(e){ #very difficult to get rid of errors here due to Error in vcov.merMod(object, use.hessian = use.hessian)... tried to write a robust function but it didn't work
					matrix(NA, nrow = 2, ncol = 4)
				})
		},
		
		compute_summary_table = function(){
			coef(summary_glm_lean(suppressWarnings(glm(private$seq_des_obj_priv_int$y ~ ., 
				data = data.frame(w = private$seq_des_obj_priv_int$w), family = "binomial"))))
		}
		
	)		
)