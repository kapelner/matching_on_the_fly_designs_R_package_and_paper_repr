#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferencePropMultiBetaRegr = R6::R6Class("SeqDesignInferencePropMultiBetaRegr",
	inherit = SeqDesignInferencePropUniBetaRegr,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 								(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 								for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),
	
	private = list(		
		generate_mod = function(){
			Xmm = private$create_design_matrix()
			# create_design_matrix is [Intercept, Treatment, Covariates]
			colnames(Xmm) = c("(Intercept)", "treatment", paste0("x", 1:(ncol(Xmm)-2)))
			
			res = fast_beta_regression_with_var(Xmm = Xmm, y = private$seq_des_obj_priv_int$y)
			
			# Ensure names are set for shared()
			names(res$b) = colnames(Xmm)
			# Standard error extraction in shared() needs vcov with names
			# We only returned ssq_b_2 (variance of treatment effect).
			# Let's return a dummy vcov matrix with the correct diagonal entry for treatment.
			p = ncol(Xmm)
			vcov_matrix = matrix(0, p, p)
			vcov_matrix[2, 2] = res$ssq_b_2
			colnames(vcov_matrix) = rownames(vcov_matrix) = colnames(Xmm)
			
			list(
				b = res$b,
				ssq_b_2 = res$ssq_b_2
			)
		}		
	)		
)
