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
		#' 								(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 								for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){		
			super$initialize(seq_des_obj, num_cores, verbose)
			assertResponseType(seq_des_obj$get_response_type(), "survival")	
		}
	),
	
	private = list(		
		generate_mod = function(){
			# Construct the full design matrix: covariates from get_X() and treatment (w)
			# The C++ function fast_weibull_regression adds its own intercept.
			# So, the matrix passed here should NOT contain an intercept.
			# private$get_X() returns the covariate matrix (without intercept) by default.
			# private$w is the treatment assignment vector.
			
			X_cov = private$get_X()
			if (ncol(X_cov) > 0) {
				full_X_matrix = cbind(X_cov, private$w)
				colnames(full_X_matrix) <- c(colnames(X_cov), "treatment")
			} else {
				full_X_matrix = matrix(private$w, ncol = 1)
				colnames(full_X_matrix) <- "treatment"
			}
			
			weibull_regr_mod = fast_weibull_regression(
				private$y,
				private$dead,
				as.matrix(full_X_matrix)
			)
			
			# Explicitly assign names to coefficients returned from C++
			# The C++ function adds an intercept, so names are (Intercept), then covariates, then treatment.
			X_names_without_intercept = colnames(full_X_matrix)
			names(weibull_regr_mod$coefficients) <- c("(Intercept)", X_names_without_intercept)

			if (is.null(weibull_regr_mod$coefficients) || is.null(weibull_regr_mod$vcov) || !is.matrix(weibull_regr_mod$vcov)){
				stop("fast_weibull_regression failed to return valid coefficients or vcov.")
			}
			
			# Construct full_coefficients and full_std_errs for consistency
			full_coefficients = c(weibull_regr_mod$coefficients, "log(scale)" = weibull_regr_mod$log_sigma)
			full_vcov = weibull_regr_mod$vcov
			colnames(full_vcov) = rownames(full_vcov) = names(full_coefficients)

			return(list(
				coefficients = full_coefficients,
				vcov = full_vcov
			))
		}
	)		
) # ADDED MISSING BRACE
