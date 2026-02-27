#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceContinMultOLS = R6::R6Class("SeqDesignInferenceContinMultOLS",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference 
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs), 
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "continuous")			
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
			full_X_matrix = private$create_design_matrix()
			mod = lm.fit(full_X_matrix, private$y) #can't beat R's built-in OLS
			private$cached_values$beta_hat_T = coef(mod)[2]
		},
		
		
		#' @description		
		#' Computes a 1-alpha level frequentist confidence interval differently for all response types, estimate types and test types.
		#' 
		#' Here we use the theory that MLE's computed for GLM's are asymptotically normal. 
		#' Hence these confidence intervals are asymptotically valid and thus approximate for any sample size.
		#' 
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' 
		#' @return 	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
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
		#' seq_des_inf$compute_mle_confidence_interval()
		#' }
		#'		
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},		
		
		
		#' @description

		
		
		#' Computes a 2-sided p-value
		#'
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' 
		#' @return 	The approximate frequentist p-value
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
		#' seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
		#' }
		#' 				
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	
	private = list(		
		compute_fast_bootstrap_distr = function(B, max_resample_attempts, n, y, dead, w) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)

			indices_mat = matrix(0L, nrow = n, ncol = B)
			
			for (b in 1:B) {
				attempt = 1
				i_b = NULL
				repeat {
					i_b = sample(n, n, replace = TRUE)
					w_b = w[i_b]
					if (any(w_b == 1, na.rm = TRUE) && any(w_b == 0, na.rm = TRUE)) {
						if (!private$any_censoring) break
						dead_b_temp = dead[i_b]
						if (any(dead_b_temp[w_b == 1] == 1) && any(dead_b_temp[w_b == 0] == 1) && min(y[i_b]) > 0) break
					}
					attempt = attempt + 1
					if (attempt > max_resample_attempts) break
				}

				if (attempt > max_resample_attempts) {
					indices_mat[, b] = rep(-1L, n) # C++ looks for -1 to return NA
				} else {
					indices_mat[, b] = i_b - 1L # Convert to 0-based for C++
				}
			}
			
			X_covars = private$get_X()
			res = compute_ols_bootstrap_parallel_cpp(y, X_covars, w, indices_mat, private$num_cores)
			return(res)
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)

			nsim = length(permutations)
			n = length(y)
			
			w_mat = matrix(0L, nrow = n, ncol = nsim)
			for (i in 1:nsim) {
				w_mat[, i] = permutations[[i]]$w
			}

			if (delta != 0) return(NULL)
			
			X_covars = private$get_X()

			res = compute_ols_distr_parallel_cpp(y, X_covars, w_mat, private$num_cores)
			return(res)
		},

		shared = function(){
#			private$cached_values$summary_table = 
#				stats::coef(summary(lm(private$seq_des_obj_priv_int$y ~ ., data = cbind(data.frame(w = private$seq_des_obj_priv_int$w), private$get_X()))))
			full_X_matrix = private$create_design_matrix()
			colnames(full_X_matrix) <- c("(Intercept)", "treatment", colnames(private$get_X()))
			
			mod = fast_ols_with_var_cpp(full_X_matrix, private$y)
			private$cached_values$beta_hat_T = mod$b[2]			
			private$cached_values$s_beta_hat_T = sqrt(mod$ssq_b_j)
			private$cached_values$is_z = FALSE 
			private$cached_values$df = private$n - (ncol(private$get_X()) + 2)
			
			# Store full coefficients and vcov for info extraction
			private$cached_values$full_coefficients = mod$b
			names(private$cached_values$full_coefficients) = colnames(full_X_matrix)
			
			# Construct vcov (only treatment var is returned by fast_ols_with_var_cpp currently, 
			# but we can at least put the treatment variance in a matrix)
			vcov_matrix = matrix(0, ncol(full_X_matrix), ncol(full_X_matrix))
			vcov_matrix[2, 2] = mod$ssq_b_j
			colnames(vcov_matrix) = rownames(vcov_matrix) = colnames(full_X_matrix)
			private$cached_values$full_vcov = vcov_matrix
		}		
	)		
)
