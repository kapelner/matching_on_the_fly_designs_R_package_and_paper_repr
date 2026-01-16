#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceIncidUnivLogRegr = R6::R6Class("SeqDesignInferenceIncidUnivLogRegr",
	inherit = SeqDesignInferenceMLEorKMforGLMs,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 								(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 								for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "incidence")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},
		
		#' @description
		#' Computes the appropriate estimate
		#' 
		#' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
		#' 
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignCRD$new(n = 6, response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceContinMultOLS$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		#' 	
		compute_treatment_estimate = function(){
			private$shared() # Ensure the model is fitted and cached values are populated
			private$cached_values$beta_hat_T
		}
	),
	
	private = list(		
		generate_mod = function(){
			# Construct the full design matrix: Intercept, covariates, and treatment
			# private$get_X() gives covariates (without intercept). private$w gives treatment.
			# Rcpp function fast_logistic_regression_cpp expects intercept, so cbind(1, ...)
			full_X_matrix = cbind(1, private$get_X(), private$w)
			
			# Add column names for clarity and easier matching with canonical models
			colnames(full_X_matrix) <- c("(Intercept)", colnames(private$get_X()), "treatment")
			
			# Call the C++ logistic regression
			logistic_regr_mod = fast_logistic_regression_cpp(full_X_matrix, private$y)
			
			# Extract coefficients
			coefficients = logistic_regr_mod$b
			names(coefficients) = colnames(full_X_matrix) # Assign names to coefficients
			
			# Compute vcov matrix (X'WX)^-1
			# W is the diagonal matrix of weights 'w' from the logistic_regr_mod
			XtWX = t(full_X_matrix) %*% diag(logistic_regr_mod$w) %*% full_X_matrix
			vcov_matrix = solve(XtWX) # Invert to get covariance matrix
			colnames(vcov_matrix) = rownames(vcov_matrix) = colnames(full_X_matrix)
			
			return(list(
				coefficients = coefficients,
				vcov = vcov_matrix
			))
		}		
	)		
)
