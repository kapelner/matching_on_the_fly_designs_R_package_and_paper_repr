#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceSurvivalUniWeibullRegr = R6::R6Class("SeqDesignInferenceSurvivalUniWeibullRegr",
	inherit = SeqDesignInferenceMLEorKMSummaryTable,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
        #' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
		#' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 							for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#' @param thin		For internal use only. Do not specify. You can thank R6's single constructor-only for this coding noise.
		#'
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE, thin = FALSE){		
			if (!thin){	
				assertResponseType(seq_des_obj$get_response_type(), "survival")			
				super$initialize(seq_des_obj, num_cores, verbose)
			}
		}
	),
	
	private = list(		
		shared = function(){
			private$compute_survival_weibull_regression(
				data_obj = private$seq_des_obj_priv_int$w
			)
		},
		
		compute_survival_weibull_regression = function(data_obj){
			tryCatch({
				surv_regr_mod = robust_survreg(private$seq_des_obj_priv_int$y, private$seq_des_obj_priv_int$dead, data_obj)	
				if (is.null(surv_regr_mod)){
					private$cached_values$summary_table = matrix(NA, nrow = 2, ncol = 4)		
				} else {
					private$cached_values$summary_table = summary(surv_regr_mod)$table
				}	
			}, error = function(e){
				private$cached_values$summary_table = matrix(NA, nrow = 2, ncol = 4)				
			})
		}		
	)		
)
		