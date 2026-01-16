#' A class that provides for relevant methods when the designs are KK matching-on-the-fly
#'
#' @description
#' An abstract class
SeqDesignInferenceKKPassThroughCompound = R6::R6Class("SeqDesignInferenceKKPassThroughCompound",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(
		
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 								(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 								for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
			private$compute_reservoir_and_match_statistics()
		}
	),
	private = list(		
		compute_reservoir_and_match_statistics = function(){	
			nRC = private$cached_values$KKstats$nRC
			nRT = private$cached_values$KKstats$nRT		
			nR = nRT + nRC #how many observations are there in the reservoir?
			m = private$cached_values$KKstats$m
					
			y_reservoir_T = private$cached_values$KKstats$y_reservoir[private$cached_values$KKstats$w_reservoir == 1] #get the reservoir responses from the treatment
			y_reservoir_C = private$cached_values$KKstats$y_reservoir[private$cached_values$KKstats$w_reservoir == 0] #get the reservoir responses from the control
			
			ssqD_bar = if (m > 1){
							var_cpp(private$cached_values$KKstats$y_matched_diffs) / m
						} else {
							NA_real_
						}
			ssqR = if (nRT > 1 && nRC > 1 && nR > 2){
						(var_cpp(y_reservoir_T) * (nRT - 1) + var_cpp(y_reservoir_C) * (nRC - 1)) / 
							(nR - 2) * (1 / nRT + 1 / nRC)
					} else {
						NA_real_
					}
			private$cached_values$KKstats$d_bar = if (m > 0) mean_cpp(private$cached_values$KKstats$y_matched_diffs) else NA_real_			
			private$cached_values$KKstats$ssqD_bar = ssqD_bar
			private$cached_values$KKstats$r_bar = if (nRT > 0 && nRC > 0) mean_cpp(y_reservoir_T) - mean_cpp(y_reservoir_C) else NA_real_
			private$cached_values$KKstats$ssqR = ssqR
			private$cached_values$KKstats$w_star = if (!is.na(ssqR) && !is.na(ssqD_bar)) ssqR / (ssqR + ssqD_bar) else NA_real_
		}
	)
)
