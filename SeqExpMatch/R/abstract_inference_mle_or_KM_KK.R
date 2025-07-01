#' Inference based on Maximum Likelihood for KK designs  
#'
#' @description
#' An abstract class
#' 
#'
SeqDesignInferenceMLEorKMKK = R6::R6Class("SeqDesignInferenceMLEorKMKK",
	inherit = SeqDesignInferenceMLEorKM,
	public = list(
		
        #' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
		#' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 							for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#'
		initialize = function(seq_des_obj, num_cores = 1, verbose = TRUE){
			if (!seq_des_obj$.__enclos_env__$private$isKK){
				stop("This type of inference is only available for KK designs.")
			}
			super$initialize(seq_des_obj, num_cores, verbose)	
			private$helper = SeqDesignInferenceHelperKK$new(seq_des_obj, private$seq_des_obj_priv_int$get_X())
		}
		
	),
	private = list(
		
	)
)
