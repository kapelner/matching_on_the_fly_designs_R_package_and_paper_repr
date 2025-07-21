#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceSurvivalUniCoxPHRegr = R6::R6Class("SeqDesignInferenceSurvivalUniCoxPHRegr",
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
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){			
			assertResponseType(seq_des_obj$get_response_type(), "survival")			
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),
	
	private = list(
		shared = function(){
			private$compute_cox_regression(
				data_obj = data.frame(w = private$seq_des_obj_priv_int$w)
			)
		},
		
		compute_cox_regression = function(data_obj){
			surv_obj = survival::Surv(private$seq_des_obj_priv_int$y, private$seq_des_obj_priv_int$dead)			
			tryCatch({
				private$cached_values$coxph_mod = suppressWarnings(survival::coxph(surv_obj ~ ., data = data_obj))
				private$cached_values$orig_summary_table = coef(summary(private$cached_values$coxph_mod))	
				private$cached_values$summary_table = matrix(NA, nrow = 2, ncol = 4)
				#convert to canonical form (dirty, but probably worth it)
				private$cached_values$summary_table[2, 1] = private$cached_values$orig_summary_table[1, 1]
				private$cached_values$summary_table[2, 2] = private$cached_values$orig_summary_table[1, 3]
				private$cached_values$summary_table[2, 4] = private$cached_values$orig_summary_table[1, 5]				
			}, error = function(e){
				private$cached_values$summary_table = matrix(NA, nrow = 2, ncol = 4)				
			})
		}		
	)		
)
